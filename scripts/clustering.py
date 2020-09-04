import random
import json
from time import time
import sys
from urllib import request
import db_utils
import argparse
import tempfile
import subprocess


def dump_lineages(db='data/gsaid.db'):
    """
    Connect to sqlite3 database and retrieve all records from the LINEAGE
    table.  Construct a dictionary keyed by accession number.
    :param db:  str, path to sqlite3 database
    :return:  dict, keyed by accession
    """
    cursor, conn = db_utils.open_connection(db)
    data = cursor.execute("SELECT `accession`, `lineage`, `probability`, `pangoLEARN_version` "
                          "FROM LINEAGE;").fetchall()
    result = {}
    for accn, lineage, prob, version in data:
        result.update({accn: {
            'lineage': lineage, 'prob': prob, 'version': version
        }})
    return result


def filter_problematic(features, url='https://raw.githubusercontent.com/W-L/ProblematicSites'
                                     '_SARS-CoV2/master/problematic_sites_sarsCov2.vcf'):
    vcf = request.urlopen(url)
    mask = {}
    for line in vcf.readlines():
        line = line.decode('utf-8')
        if line.startswith('#'):
            continue
        _, pos, _, ref, alt, _, filt, info = line.strip().split()
        if filt == 'mask':
            mask.update({int(pos)-1: {  # convert to 0-index
                'ref': ref, 'alt': alt, 'info': info}
            })

    # apply filters to feature vectors
    count = 0
    for row in features:
        filtered = []
        for typ, pos, alt in row['diffs']:
            if typ == '~' and int(pos) in mask and alt in mask[pos]['alt']:
                continue
            if typ != '-' and 'N' in alt:
                # drop substitutions and insertions with uncalled bases
                continue
            filtered.append(tuple([typ, pos, alt]))

        count += len(row['diffs']) - len(filtered)
        row['diffs'] = filtered

    print('filtered {} problematic features'.format(count))
    return features


def sample_with_replacement(population, k):
    n = len(population)
    pop = list(population)
    return [pop[int(n * random.random())] for _ in range(k)]


def manhattan_dist(fv1, fv2, boot=None):
    """
    Calculate the Manhattan distance between two feature vectors
    (differences of genomes from reference).

    :param fv1:  list, feature vector for first genome (differences
                 relative to reference)
    :param fv2:  list, feature vector for second genome
    :param un:  set, union of features (bootstrapping)
    """
    symdiff = set(fv1) ^ set(fv2)  # symmetric difference
    if boot:
        return sum([boot.count(feat) for feat in symdiff])
    else:
        return len(symdiff)


def total_missing(row):
    res = 0
    for left, right in row['missing']:
        res += right-left
    return res


def import_json(path, max_missing=600):
    with open(path) as fp:
        features = json.load(fp)

    # remove genomes with too many uncalled bases
    count = len(features)
    features = [row for row in features if total_missing(row) < max_missing]
    print("dropped {} records with uncalled bases in excess of {}".format(
        count - len(features), max_missing
    ))

    # remove features known to be problematic
    features = filter_problematic(features)

    return features


def apply_features(row, refseq):
    """
    Reconstitute genome sequence from feature vector (genetic differences) and
    missing data vector.
    """
    pass


def split_by_lineage(features, lineages):
    """
    Partition feature list by Pangolin lineage assignments
    :param features:  list, feature vectors by GISAID record
    :param lineages:  dict, lineage assignment keyed by accession, from dump_lineages()
    :return:  dict, feature lists keyed by lineage
    """
    result = {}
    for row in features:
        accn = row['name'].split('|')[1]
        val = lineages.get(accn, None)
        if val is None:
            print("Error in clustering::split_by_lineage(), no lineage assignment"
                  " for accession {}".format(accn))
            sys.exit()
        lineage = val['lineage']

        if lineage not in result:
            result.update({lineage: []})
        result[lineage].append(row)
    return result


def get_union(features):
    un = set()
    for row in features:
        un = un.union(set(row['diffs']))
    return un


def write_dists(outfile, features):
    """
    Write pairwise distance matrix in PHYLIP format to file.
    If list contains fewer than 3 feature vectors, then return False.
    :param outfile:  file, open stream in write mode
    :param features:  list, feature vectors
    :return:  True if successful
    """
    n = len(features)
    if n < 3:
        # not much point making this tree
        return False

    # print('generating union')
    # un = get_union(features)
    # number of taxa
    outfile.write('{0:>5}\n'.format(n))
    for i in range(n):
        row = features[i]
        outfile.write('{}'.format(row['name']))  # .split('_')[-1]))

        for j in range(n):
            md = manhattan_dist(features[i]['diffs'], features[j]['diffs'])
            outfile.write(' {0:>2}'.format(md))
        outfile.write('\n')
    return True


def rapidnj(features, negative=False, binpath='rapidnj'):
    """
    Wrapper function for rapidNJ.  Writes a pairwise distance matrix to a
    temporary file and calls rapidNJ to apply neighbor-joining to the matrix.

    :param features:  list, differences from reference as feature vectors
    :param negative:  if True, allow negative branch lengths
    :return:  str, Newick tree string
    """
    with tempfile.NamedTemporaryFile('w', delete=False) as outfile:
        success = write_dists(outfile, features)

    if not success:
        return None

    cmd = [binpath, outfile.name, '-i', 'pd']
    if not negative:
        cmd.append('--no-negative-length')
    stdout = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
    return stdout.decode('utf-8')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Partition genomes by Pangolin lineage, generate distance"
                    "matrices and use neighbor-joining."
    )
    parser.add_argument("json", type=str,
                        help="input, JSON file generated by minimap2.py")
    parser.add_argument("--db", type=str, default="data/gsaid.db",
                        help="Path to sqlite3 database")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    print ('loading lineage classifications from database')
    lineages = dump_lineages(args.db)

    print('loading JSON')
    features = import_json(args.json)

    by_lineage = split_by_lineage(features, lineages)
    for lineage, lfeatures in by_lineage.items():
        nwk_str = rapidnj(lfeatures)
        if nwk_str is None:
            continue

        with open('{}.nwk'.format(lineage), 'w') as outfile:
            outfile.write(nwk_str)

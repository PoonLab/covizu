import random
import json
from time import time
import sys
from urllib import request
import db_utils
import argparse
import tempfile
import subprocess
from functools import reduce
from Bio import Phylo
from io import StringIO


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


def total_missing(row):
    res = 0
    for left, right in row['missing']:
        res += right-left
    return res


def import_json(path, max_missing=600):
    with open(path) as fp:
        features = json.load(fp)

    # remove features known to be problematic
    features = filter_problematic(features)

    # remove genomes with too many uncalled bases
    count = len(features)
    features = [row for row in features if total_missing(row) < max_missing]
    print("dropped {} records with uncalled bases in excess of {}".format(
        count - len(features), max_missing
    ))

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


def sample_with_replacement(population, k):
    """ Sample <k> items from iterable <population> with replacement """
    n = len(population)
    pop = list(population)
    return [pop[int(n * random.random())] for _ in range(k)]


def get_distances(features, nboot=1):
    """
    Write pairwise distance matrix in PHYLIP format to file.

    :param features:  list, feature vectors from load_json()/split_by_lineage()
    :param nboot:  int, number of bootstrap samples
    :return:  list, nested lists representing <n x n x nboot> matrix of
              pairwise distances
    """
    n = len(features)

    # recast feature vectors as sets
    fvs = [set(f['diffs']) for f in features]

    # calculate symmetric differences
    sym_diffs = {}
    for i in range(n-1):
        for j in range(i+1, n):
            sym_diffs[(i, j)] = fvs[i] ^ fvs[j]

    # bootstrap sampling from feature set union
    dists = [[[0]*nboot for _ in range(n)] for _ in range(n)]
    union = reduce(lambda x, y: x.union(y), fvs)
    for k in range(nboot):
        sample = sample_with_replacement(union, len(union))
        weights = dict([(f, sample.count(f)) for f in set(sample)])
        for i in range(n-1):
            for j in range(i+1, n):
                dists[j][i][k] = dists[i][j][k] = sum([
                    weights.get(f, 0) for f in sym_diffs[(i, j)]
                ])
    return dists


def write_dists(k, dists, labels):
    """
    Export pairwise distance matrix to a temporary file.
    :param k:  int, bootstrap sample number
    :param dists:  list, nested lists of pairwise distances from get_distances()
    :param features:  list, sequence names and feature vectors
    :return:  file object, temporary file
    """
    outfile = tempfile.NamedTemporaryFile('w', delete=False)
    n = len(dists)

    outfile.write('{0:>5}\n'.format(n))
    for i in range(n):
        outfile.write('{}'.format(features[i]['name']))
        for j in range(n):
            outfile.write(' {0:>2}'.format(dists[i][j][k]))
        outfile.write('\n')

    outfile.close()
    return outfile


def rapidnj(dists, labels, negative=False, binpath='rapidnj'):
    """
    Wrapper function for rapidNJ.  Writes a pairwise distance matrix to a
    temporary file and calls rapidNJ to apply neighbor-joining to the matrix.

    :param features:  list, differences from reference as feature vectors
    :param negative:  if True, allow negative branch lengths
    :yield:  Biopython.Phylo.BaseTree objects
    """
    nboot = len(dists[0][0])
    for k in range(nboot):
        outfile = write_dists(k, dists, labels)

        cmd = [binpath, outfile.name, '-i', 'pd']
        if not negative:
            cmd.append('--no-negative-length')

        stdout = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
        with StringIO(stdout.decode('utf-8')) as handle:
            yield Phylo.read(handle, 'newick')


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

    print('loading lineage classifications from database')
    lineages = dump_lineages(args.db)

    print('loading JSON')
    features = import_json(args.json)

    by_lineage = split_by_lineage(features, lineages)
    for lineage, lfeatures in by_lineage.items():
        print(lineage)
        dists = get_distances(lfeatures, 10)
        labels = [row['name'] for row in lfeatures]

        result = [phy for phy in rapidnj(dists, labels)]
        for i, phy in enumerate(result):
            Phylo.write(phy, file='{}-{}.nwk'.format(lineage, i), format='newick')
        break

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
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.

    :param features:  list, return object from import_json()
    :param url:  str, URL to the VCF file
    :return:
    """
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
    """ Calculate the total number of missing sites from closed-open interval annotations """
    res = 0
    for left, right in row['missing']:
        res += right-left
    return res


def import_json(path, max_missing=600):
    """
    Read genome features (genetic differences from reference) from JSON file.

    :param path:  str, relative or absolute path to JSON input file
    :param max_missing:  int, maximum tolerated number of uncalled bases
    :return:  list, each entry comprises a dict with keys 'name', 'diffs', and 'missing'
    """
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

    :param row:  dict, entry from features list returned by import_json()
    :param refseq:  str, reference genome
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

    :param dists:  list, an <n x n x B> matrix as nested lists, where `n` is the
                   number of genomes and `B` is the number of bootstrap samples.
    :param labels:  list, genome names for labeling tree tips.
    :param negative:  if True, allow negative branch lengths
    :param binpath:  str, path to rapidNJ executable

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


def consensus_tree(iter, minboot=0.5):
    """
    Generate consensus tree from trees reconstructed from bootstrap samples.
    FIXME: Bio.Phylo.Consensus has function majority_consensus() - is this for rooted trees only?

    :param iter:  generator, Phylo.BaseTree objects from rapidnj()
    :param minboot:  float, minimum threshold to maintain clade/split
    :return:
    """
    splits = {}
    for phy in iter:
        tips = [tip.name for tip in phy.get_terminals()]
        dups = [tip for tip in set(tips) if tips.count(tip) > 1]
        if dups:
            print("ERROR: in consensus_tree(), found duplicate tip labels:")
            for tipname in dups:
                print(tipname)
            sys.exit()

        all_tips = set(tips)
        for clade in phy.get_nonterminals():
            if clade is phy.root:
                continue
            my_tips = set([tip.name for tip in clade.get_terminals()])
            their_tips = all_tips - my_tips

            key = (my_tips, their_tips) if len(my_tips) < len(their_tips) else (their_tips, my_tips)
            if key not in splits:
                splits.update({key: 0})
            splits[key] += 1

    # construct consensus tree from majority splits




def parse_args():
    """ Command-line interface """
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

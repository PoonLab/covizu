import networkx as nx
from datetime import date
import sys
import argparse
from gotoh2 import iter_fasta
import csv

MIN_DIST_CUTOFF = 0.0005  # ~15 nt differences from any other genome


def parse_label(label):
    """
    Extract country and date of sample collection from GISAID label.
    Return date as None if collection date is ambiguous.

    :param label:  str, GISAID sequence label
    :return: (country, date)
    """
    info, epi_id, ymd = label.split('|')
    country = info.split('/')[1]
    try:
        year, month, day = list(map(int, ymd.split('-')))
        return country, date(year, month, day)
    except ValueError:
        return country, None
    except:
        raise


def check_cliques(node, cliques):
    for idx, clique in enumerate(cliques):
        if node in clique:
            return idx
    return None


def clustering(tn93_file, callback=None):
    """
    Use TN93 distances:
     tn93 -o data/gisaid.tn93.csv data/gisaid-filtered.fa

    to cluster genome sequences into unique variants.  Simultaneously,
    exclude sequences that are too distant from any other sequence.

    @param tn93_file: input, path to CSV file with TN93 distances
    @param callback: optional, function to pass messages to stdout
    @return networkx Graph object
    """

    handle = open(tn93_file)
    _ = next(handle)

    if callback:
        callback("building graph from nodes to find clusters")

    # filter edges that are below threshold distance
    stars = {}
    nodes = set()
    for line in handle:
        id1, id2, dist = line.strip().split(',')
        dist = float(dist)
        for node in [id1, id2]:
            if node not in nodes:
                nodes.update({node})

        if dist < 1e-09:
            if id1 not in stars:
                stars.update({id1: set()})
            stars[id1].update({id2})
            if id2 not in stars:
                stars.update({id2: set()})
            stars[id2].update({id1})

    handle.close()

    # assemble cliques from stars - takes a few minutes
    cliques = []
    for ego, alters in stars.items():
        clique = alters.union({ego})
        for alter in alters:
            # alters of alter must be entirely contained within clique, minus self
            if clique.difference(stars[alter]) != {alter}:
                clique.remove(alter)
        cliques.append(clique)

    # generate unique set of cliques (frozensets are hashable)
    uniques = list(set([frozenset(c) for c in cliques]))

    # generate the intersection graph
    ig = nx.Graph()
    for i in range(len(uniques)):
        cliq1 = uniques[i]
        # node is weighted by number of sequences
        ig.add_node(i, weight=len(cliq1))
        for j in range(i, len(uniques)):
            cliq2 = uniques[j]
            if len(cliq1.intersection(cliq2)) == 0:
                ig.add_edge(i, j)

    # find maximal clique in intersection graph with highest total weight
    # FIXME: this is the most demanding step - need to explore parallel options
    max_weight = 0
    max_sg = None
    for igc in nx.find_cliques(ig):
        sg = ig.subgraph(igc)
        weights = nx.get_node_attributes(sg, 'weight').values()
        total_weight = sum(list(weights))
        if total_weight > max_weight:
            max_weight = total_weight
            max_sg = sg
            print(max_weight)

    return max_sg



def write_variants(G, csv_file, fasta_in, fasta_out, callback=None):
    """
    Write CSV file describing the content of each genome variant cluster.
    :param G:  networkx.graph object from clustering()
    :param csv_file:  path to write variants in CSV format
    :param fasta_in:  path to FASTA file of aligned genomes
    :param fasta_out:  path to write FASTA file of unique genome variants
    :param callback:  optional, for passing messages to stdout
    :return: dict, {label: collection date}
    """
    components = list(nx.connected_components(G))
    if callback:
        callback("graph comprises {} sequences and "
          "{} components".format(len(G), len(components)))

    # generate cluster information
    writer = csv.writer(open(csv_file, 'w'))
    writer.writerow(['cluster', 'label', 'coldate', 'country'])

    clusters = {}
    for cluster in components:
        subG = G.subgraph(cluster)
        # omit records with ambiguous collection dates
        intermed = [
            (ndata['coldate'], ndata['country'], node) for node, ndata
            in subG.nodes(data=True) if ndata['coldate'] is not None
        ]
        intermed.sort()  # increasing order of collection dates
        if len(intermed) == 0:
            # exclude clusters with all missing dates
            continue

        # label cluster by earliest case
        _, _, label = intermed[0]
        clusters.update({label: len(cluster)})

        # write cluster contents to info file
        for coldate, country, node in intermed:
            writer.writerow([label, node, coldate, country])

    outfile = open(fasta_out, 'w')
    for h, s in iter_fasta(open(fasta_in)):
        h = h.strip()
        if h not in clusters:
            continue
        outfile.write(">{}\n{}\n".format(h, s.replace('?', 'N')))

    return clusters


def parse_args():
    parser = argparse.ArgumentParser(
        description="Processing and clustering of aligned SARS-CoV-2"
                    " genome sequences into unique variants."
    )
    parser.add_argument('--tn93', default='data/gisaid.tn93.csv',
                        help='input, path to CSV file containing TN93 '
                             'distances.')
    parser.add_argument('--csv_out', default='data/variants.csv',
                        help='output, path to write CSV describing '
                             'composition of variants')
    parser.add_argument('--fasta_in', default='data/gisaid-filtered.fa',
                        help='input, path to FASTA with aligned genomes')
    parser.add_argument('--fasta_out', default='data/variants.fa',
                        help='output, path to write cluster FASTA')

    return parser.parse_args()


if __name__ == "__main__":
    def callback(msg):
        sys.stdout.write(msg+'\n')
        sys.stdout.flush()

    args = parse_args()
    G = clustering(args.tn93, callback=callback)
    write_variants(G, csv_file=args.csv_out,
               fasta_in=args.fasta_in, fasta_out=args.fasta_out,
               callback=callback)

import networkx as nx
from networkx.algorithms import tree
import argparse
from csv import DictReader
import sys
from gotoh2 import iter_fasta


def build_mst(tn93_file, callback=None):
    """
    :param tn93_file:  path to CSV of TN93 distances
    :param callback:  optional, progress monitoring function
    """
    if callback:
        callback("populate complete graph with TN93 distances")
    G = nx.Graph()
    with open(tn93_file) as f:
        _ = next(f)
        for line in f:
            id1, id2, dist = line.strip().split(',')
            dist = float(dist)
            for node in [id1, id2]:
                if node not in G:
                    G.add_node(node)
            G.add_edge(id1, id2, weight=dist)

    if callback:
        callback("generating minimum spanning tree")
    return tree.minimum_spanning_tree(G, weight='weight')


def get_edgelist(g):
    """
    Generate dictionary of edge weights for rapid access.
    :param g:  networkx.Graph object
    :return:  dict
    """
    edgelist = {}
    for n1, n2, dist in g.edges(data='weight'):
        if n1 not in edgelist:
            edgelist.update({n1: {}})
        edgelist[n1].update({n2: dist})
        if n2 not in edgelist:
            edgelist.update({n2: {}})
        edgelist[n2].update({n1: dist})
    return edgelist


def parse_clusters(infile):
    """
    :param infile: cluster info file from clustering.py
    :return: dict
    """
    clusters = {}
    reader = DictReader(open(infile))
    for row in reader:
        cluster = row['cluster']
        if cluster not in clusters:
            clusters.update({cluster: {
                'count': 0, 'coldate': row['coldate']
            }})
        clusters[cluster]['count'] += 1
    return clusters


def traversal(node, parent, edgelist, history):
    """
    Traverse minimum spanning tree from earliest variant.
    Search edge list for children of root node and recurse.
    :param node:  current node in recursion
    :param parent:  upstream node
    :param edgelist:  dict, from get_edgelist()
    :param history:  list, track nodes visited in graph
    :yield:  nodes by preorder traversal
    """
    history.append(node)
    yield node, parent
    children = [child for child, _ in edgelist[node].items()
                if child not in history]
    for child in children:
        for obj in traversal(child, node, edgelist, history=history):
            yield obj


def cut_tree(mst, clusters, outstem, cutoff=15):
    """
    Root the minimum spanning tree on the earliest node and
    traverse the tree.  Remove in-edges for nodes with out-degree
    exceeding cutoff.
    :param mst:  networkx.Graph
    :param clusters:  dict, object returned from parse_clusters()
    :param outstem:  formatted string for writing output files with
                     a single integer index
    :param cutoff:  out-degree threshold for cutting subtrees; because
                    input is a tree, in-degree is always 1.
    """
    subtrees = []

    # root the MST on the cluster with most cases
    intermed = [(nd['count'], node) for node, nd in clusters.items()]
    intermed.sort(reverse=True)
    root = intermed[0][1]

    edgelist = get_edgelist(mst)
    dg = nx.DiGraph()
    for child, parent in traversal(root, None, edgelist, history=[]):
        if parent is None:
            continue
        if mst.degree[child] > cutoff:
            # omit edge to initialize new subtree
            subtrees.append(child)
            continue
        dg.add_edge(parent, child)

    # export components
    components = list(nx.weakly_connected_components(dg))
    for i, comp in enumerate(components):
        with open(outstem.format(i), 'w') as outfile:
            outfile.write('parent,child,dist\n')
            # extract subtree
            sg = nx.subgraph(dg, comp)
            for parent, child in sg.edges():
                dist = edgelist[parent][child]
                outfile.write('{},{},{}\n'.format(parent, child, dist))

    return subtrees


def write_treetime_inputs(clusters, subtrees, n_output, fasta_in,
                          fasta_out, dates_out):
    """
    Export cluster sequences for TimeTree analysis
    :param clusters:  dict, returned from write_info()
    :param subtrees:  list, labels for nodes that root subtrees in MST
    :param n_output:  int, number of cluster sequences to output
    :param fasta_in:  input, path to FASTA containing genome sequences
    :param fasta_out: output, path to write FASTA of cluster sequences
    :param dates_out: output, path to write cluster dates for TreeTime
    """
    # filter the largest clusters
    intermed = [(v['count'], k) for k, v in clusters.items()]
    intermed.sort(reverse=True)
    keys = [key for count, key in intermed[:n_output]]
    keys += [node for node in subtrees if node not in keys]

    # open file streams
    outfile = open(fasta_out, 'w')
    datefile = open(dates_out, 'w')
    datefile.write('name,date\n')

    for h, s in iter_fasta(open(fasta_in)):
        if h not in keys:
            continue
        accession = h.split('|')[1]
        outfile.write(">{}\n{}\n".format(accession, s.replace('?', 'N')))
        datefile.write('{},{}\n'.format(accession, clusters[h]['coldate']))

    outfile.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tn93', default='data/clusters.tn93.csv',
                        help='input, path to TN93 CSV file')
    parser.add_argument('--info', default='data/clusters.info.csv',
                        help='input, path to CSV with cluster information, '
                             'generated by clustering.py')
    parser.add_argument('--outstem', default='mst/component-{}.edgelist.csv',
                        help='output, stem for output files with Python '
                             'formatted string syntax with one placeholder '
                             '("{}") for cluster index.')
    parser.add_argument('--cutoff', type=int, default=15,
                        help='option, threshold out-degree to break a subtree'
                             ' from the minimum spanning tree; default 15')

    parser.add_argument('--infasta', default='data/clusters.fa',
                        help='Path to FASTA file with aligned GISAID data.')
    parser.add_argument('--outfasta', default='data/treetime.fa',
                        help='Path to file to write FASTA output.')
    parser.add_argument('--dates', default='data/treetime.csv',
                        help='Path to file to write dates in CSV format '
                             '(for TreeTime analysis).')
    parser.add_argument('-n', default=20,
                        help='Number of sequences to export to TreeTime.')

    return parser.parse_args()


if __name__ == '__main__':
    def callback(msg):
        sys.stdout.write(msg+'\n')
        sys.stdout.flush()

    args = parse_args()
    clusters = parse_clusters(args.info)
    mst = build_mst(args.tn93, callback=callback)
    subtrees = cut_tree(mst, clusters, outstem=args.outstem,
                        cutoff=args.cutoff)
    write_treetime_inputs(clusters, subtrees, n_output=args.n,
                          fasta_in=args.infasta, fasta_out=args.outfasta,
                          dates_out=args.dates)

# tn93 -t 0.0002 -o data/clusters.tn93.csv data/clusters.fa

# pass outputs to fasttree2 and treetime
# fasttree2 -nt < clusters.fa > clusters.ft2.nwk
# python3 prune-long-tips.py
# treetime --tree data/clusters.pruned.nwk --aln data/clusters.fa --dates data/clusters.dates.csv
# python3 parse-nexus.py
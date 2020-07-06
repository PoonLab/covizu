import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
from datetime import date
import sys
import argparse
from gotoh2 import iter_fasta
import csv


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


def import_graph(tn93_file, mindist=1e-09, callback=None):
    """
    Use TN93 distances:
      tn93 -o data/gisaid.tn93.csv data/gisaid-filtered.fa
    to generate a graph where an edge connects sequences with a
    distance of effectively zero ("identical").

    @param tn93_file: input, path to CSV file with TN93 distances
    @param callback: optional, function to pass messages to stdout
    @return networkx Graph object
    """
    if callback:
        callback("importing TN93 distances")
    graph = nx.Graph()
    with open(tn93_file) as handle:
        _ = next(handle)  # skip header line
        for line in handle:
            id1, id2, dist = line.strip().split(',')
            for node in [id1, id2]:
                if node not in graph:
                    graph.add_node(node)

            if float(dist) < mindist:
                graph.add_edge(id1, id2)
    if callback:
        callback("built graph with {} nodes".format(len(graph)))
    return graph


def clique_clustering(graph):
    """
    Use maximal cliques to define variants, where a
    maximal clique is the largest subgraph containing a given node
    such that every pair of nodes is connected by an edge.
    WARNING: this is a very time-consuming process, not recommended
    for large graphs!

    :param graph: networkx.Graph object from import_graph()
    :return:  list of connected components as lists of node labels
    """
    result = []
    for component in nx.connected_components(graph):
        sg = graph.subgraph(component)
        cliques = nx.find_cliques(sg)

        # generate unique set of cliques (frozensets are hashable)
        uniques = list(set([frozenset(c) for c in cliques]))

        # generate clique graph where edges indicate non-overlapping cliques
        cgraph = nx.Graph()
        for i in range(len(uniques)):
            cliq1 = uniques[i]
            # node is weighted by number of sequences in clique
            cgraph.add_node(i, weight=len(cliq1))
            for j in range(i, len(uniques)):
                cliq2 = uniques[j]
                if len(cliq1.intersection(cliq2)) == 0:
                    if j not in cgraph:
                        cgraph.add_node(j, weight=len(cliq2))
                    cgraph.add_edge(i, j)

        # find maximal clique in clique graph with highest total weight
        max_weight = 0
        max_csg = None
        for cclique in nx.find_cliques(cgraph):
            # FIXME: some iterations may be redundant
            csg = cgraph.subgraph(cclique)
            weights = nx.get_node_attributes(csg, 'weight').values()
            total_weight = sum(list(weights))
            if total_weight > max_weight:
                max_weight = total_weight
                max_csg = csg

        # generate partition of node labels
        for clique in max_csg.nodes():
            clique_set = [uniques[i] for i in clique]
            result.extend(clique_set)

    return result


def modularity_clustering(graph, size_cutoff=10, deg_cutoff=0.5, 
                          callback=None):
    """
    Use the Clauset-Newman-Moore greedy modularity maximization
    algorithm to partition the TN93 pairwise graph into communities.
    Modularity quantifies the density of edges at the periphery of
    a community relative to the density within it.
    TODO: try other methods like Louvain algorithm

    :param graph:  networkx.Graph object from import_graph()
    :param size_cutoff:  int, minimum component size to consider
                         applying modularity community detection
    :param deg_cutoff:  float, maximum edge density at which use
                        community detection.
    :param callback:  optional, write verbose messages
    :return: list, lists of node labels
    """
    result = []
    for component in nx.connected_components(graph):
        if len(component) > size_cutoff:
            sg = graph.subgraph(component)
            # retrieve list of degree sizes
            deg = [d for _, d in sg.degree()]
            mean_deg = sum(deg) / float(len(deg))
            if mean_deg / len(deg) < deg_cutoff:
                communities = list(greedy_modularity_communities(sg))
                if callback:
                    callback(
                        'Partitioning component of size {} into {}'
                        'communities'.format(len(component), len(communities))
                    )
                result.extend(communities)
            else:
                # component has sufficient edge density
                result.append(component)
        else:
            result.append(component)

    return result


def write_variants(components, csv_file, fasta_in, fasta_out, callback=None):
    """
    Write CSV file describing the content of each genome variant cluster.
    :param components:  list, lists of node labels
    :param csv_file:  str, path to write variants in CSV format
    :param fasta_in:  str, path to FASTA file of aligned genomes
    :param fasta_out:  str, path to write FASTA file of unique genome variants
    :param callback:  optional, for passing messages to stdout
    :return: dict, {label: collection date}
    """

    # generate cluster information
    writer = csv.writer(open(csv_file, 'w'))
    writer.writerow(['cluster', 'label', 'coldate', 'country'])

    variants = {}
    for component in components:
        # omit records with ambiguous collection dates
        intermed = []
        for label in component:
            country, coldate = parse_label(label)
            if coldate is None:
                continue
            intermed.append([coldate, country, label])

        if len(intermed) == 0:
            # exclude components with all missing dates
            continue

        intermed.sort()  # increasing order of collection dates

        # label cluster by earliest case
        _, _, label0 = intermed[0]
        variants.update({label0: len(component)})

        # write cluster contents to info file
        for coldate, country, label in intermed:
            writer.writerow([label0, label, coldate, country])

    # write reduced FASTA file
    outfile = open(fasta_out, 'w')
    for h, s in iter_fasta(open(fasta_in)):
        h = h.strip()
        if h not in variants:
            continue
        outfile.write(">{}\n{}\n".format(h, s.replace('?', 'N')))

    return variants


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
    graph = import_graph(args.tn93, callback=callback)
    components = modularity_clustering(graph, callback=callback)
    write_variants(components, csv_file=args.csv_out,
                   fasta_in=args.fasta_in, fasta_out=args.fasta_out,
                   callback=callback)

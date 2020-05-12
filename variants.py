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

def read_countries(country_file):
    """
    :param country_file: path to CSV of country, continent
    :return dict
    """
    # load country to region mapping
    countries = {}
    with open(country_file) as handle:
        for line in handle:
            try:
                country, region = line.strip().split(',')
            except:
                print(line)
                raise
            countries.update({country: region})
    return countries


def clustering(tn93_file, countries, callback=None):
    """
    Use TN93 distances:
     tn93 -o data/gisaid.tn93.csv data/gisaid-filtered.fa

    to cluster genome sequences into unique variants.  Simultaneously,
    exclude sequences that are too distant from any other sequence.

    @param tn93_file: input, path to CSV file with TN93 distances
    @param countries: dict, from read_countries()
    """
    if callback:
        callback("Filter sequences that are too distant")
    min_dists = {}
    unknown_country = []
    handle = open(tn93_file)
    _ = next(handle)
    for line in handle:
        id1, id2, dist = line.strip().split(',')
        dist = float(dist)
        for node in [id1, id2]:
            if node not in min_dists:
                min_dists.update({node: 10.0})
                # check for new country values while we're at it
                country, coldate = parse_label(node)
                if country not in countries:
                    unknown_country.append(country)
            if dist < min_dists[node]:
                min_dists[node] = dist

    if unknown_country:
        print('Error: unrecognized <country> annotations detected:')
        for country in set(unknown_country):
            print(country)
        print("Please update the countries.csv file.")
        sys.exit()

    intermed = [(dist, node) for node, dist in min_dists.items()]
    intermed.sort(reverse=True)
    discard = {}
    for dist, node in intermed:
        if dist < MIN_DIST_CUTOFF:
            # since the list is sorted, we can stop search
            break
        discard.update({node: dist})

    if discard and callback:
        callback("flagged {} sequences for discarding".format(len(discard)))
        for node, dist in discard.items():
            callback('{} ({})'.format(node, dist))

    if callback:
        callback("building graph from nodes to find clusters")
    handle.seek(0)  # reset TN93 file
    _ = next(handle)
    G = nx.Graph()
    for line in handle:
        id1, id2, dist = line.strip().split(',')
        dist = float(dist)
        if id1 in discard or id2 in discard:
            continue
        for node in [id1, id2]:
            if node not in G:
                country, coldate = parse_label(node)
                G.add_node(node, country=country, coldate=coldate)

        if dist < 1e-09:
            # some very small distances reported - mixtures?
            G.add_edge(id1, id2)

    return G


def write_info(G, countries, info_file, fasta_in, fasta_out, callback=None):
    """
    Write CSV file describing the content of each genome variant cluster.
    :param G:  networkx.graph object from clustering()
    :param countries:  dict, mapping country to region (continent)
    :param info_file:  output file path
    :return: dict, label: collection date
    """
    components = list(nx.connected_components(G))
    if callback:
        callback("graph comprises {} sequences and "
          "{} components".format(len(G), len(components)))

    # generate cluster information
    writer = csv.writer(open(info_file, 'w'))
    writer.writerow(['cluster', 'label', 'coldate', 'region', 'country'])

    clusters = {}
    for cluster in components:
        subG = G.subgraph(cluster)
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
            region = countries[country]
            writer.writerow([label, node, coldate, region, country])

    outfile = open(fasta_out, 'w')
    for h, s in iter_fasta(open(fasta_in)):
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
    parser.add_argument('--country', default='countries.csv',
                        help='input, path to CSV file linking countries '
                             'to geographic regions (continents).')
    parser.add_argument('--info', default='data/variants.csv',
                        help='output, path to write CSV describing '
                             'composition of variants')
    parser.add_argument('--fasta_in', default='data/gisaid-aligned.fa',
                        help='input, path to FASTA with aligned genomes')
    parser.add_argument('--fasta_out', default='data/variants.fa',
                        help='output, path to write cluster FASTA')

    return parser.parse_args()


if __name__ == "__main__":
    def callback(msg):
        sys.stdout.write(msg+'\n')
        sys.stdout.flush()

    args = parse_args()

    countries = read_countries(args.country)
    G = clustering(args.tn93, countries=countries, callback=callback)
    write_info(G, countries=countries, info_file=args.info,
               fasta_in=args.fasta_in, fasta_out=args.fasta_out,
               callback=callback)

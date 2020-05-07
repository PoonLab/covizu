import networkx as nx
from datetime import date
import sys
import argparse
from gotoh2 import iter_fasta

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


def clustering(tn93_file, country_file, callback=None):
    """
    Use TN93 distances:
     tn93 -o data/gisaid.tn93.csv data/gisaid-filtered.fa

    to cluster genome sequences into unique variants.  Simultaneously,
    exclude sequences that are too distant from any other sequence.

    @param country_file: input, path to CSV file with country data
    @param tn93_file: input, path to CSV file with TN93 distances
    """

    # load country to region mapping
    country2region = {}
    with open(country_file) as handle:
        for line in handle:
            try:
                country, region = line.strip().split(',')
            except:
                print(line)
                raise
            country2region.update({country: region})

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
                if country not in country2region:
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

    return G, country2region


def write_info(G, country2region, info_file, callback=None):
    """
    Write CSV file describing the content of each genome variant cluster.
    :param G:  networkx.graph object from clustering()
    :param country2region:  dict, mapping country to region (continent)
    :param info_file:  output file path
    :return: dict, label: collection date
    """
    components = list(nx.connected_components(G))
    if callback:
        callback("graph comprises {} sequences and "
          "{} components".format(len(G), len(components)))

    # generate cluster information
    infofile = open(info_file, 'w')
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
        coldate, _, label = intermed[0]
        clusters.update({label: {'coldate': coldate, 'count': len(cluster)}})

        # write cluster contents to info file
        for coldate, country, node in intermed:
            region = country2region[country]
            infofile.write('{},{},{},{},{}\n'.format(
                label, node, coldate, region, country
            ))
    return clusters


def write_treetime_inputs(k, clusters, fasta_in, fasta_out, dates_out):
    """
    Export cluster sequences for TimeTree analysis
    :param k:  int, number of cluster sequences to output
    :param clusters:  dict, returned from write_info()
    :param fasta_in:  input, path to FASTA containing genome sequences
    :param fasta_out: output, path to write FASTA of cluster sequences
    :param dates_out: output, path to write cluster dates for TreeTime
    """

    # filter the largest clusters
    intermed = [(v['count'], k) for k, v in clusters.items()]
    intermed.sort(reverse=True)
    keys = [key for count, key in intermed[:k]]

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


# pass outputs to fasttree2 and treetime
# fasttree2 -nt < clusters.fa > clusters.ft2.nwk
# python3 prune-long-tips.py
# treetime --tree data/clusters.pruned.nwk --aln data/clusters.fa --dates data/clusters.dates.csv
# python3 parse-nexus.py


def parse_args():
    parser = argparse.ArgumentParser(
        description="Processing and clustering of aligned SARS-CoV-2"
                    " genome sequences into unique variants."
    )
    parser.add_argument('--country', default='countries.csv',
                        help='Path to CSV file linking countries to '
                             'regions.')
    parser.add_argument('--tn93', default='data/gisaid.tn93.csv',
                        help='Path to CSV file containing TN93 '
                             'distances.')
    parser.add_argument('--info', default='data/clusters.info.csv',
                        help='Path to CSV containing cluster info '
                             '(generated by clustering.py)')
    parser.add_argument('--infasta', default='data/gisaid-aligned.fa',
                        help='Path to FASTA file with aligned GISAID data.')
    parser.add_argument('--outfasta', default='data/clusters.fa',
                        help='Path to file to write FASTA output.')
    parser.add_argument('--dates', default='data/clusters.dates.csv',
                        help='Path to file to write dates in CSV format '
                             '(for TreeTime analysis).')
    parser.add_argument('-k', default=20,
                        help='Number of sequences to export to TreeTime.')

    return parser.parse_args()


if __name__ == "__main__":
    def callback(msg):
        sys.stdout.write(msg)
        sys.stdout.flush()

    args = parse_args()
    G, c2r = clustering(
        tn93_file=args.tn93, country_file=args.country, callback=callback
    )
    clusters = write_info(G, country2region=c2r, info_file=args.info,
                          callback=callback)
    write_treetime_inputs(clusters, fasta_in=args.infasta, fasta_out=args.outfasta,
                          dates_out=args.dates)

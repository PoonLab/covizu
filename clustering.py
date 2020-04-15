import networkx as nx
from networkx.algorithms import tree
from datetime import date
from csv import DictWriter
import math
import sys
from gotoh2 import iter_fasta

MIN_DIST_CUTOFF = 0.0005

def parse_label(label):
    info, epi_id, ymd = label.split('|')
    country = info.split('/')[1]
    try:
        year, month, day = list(map(int, ymd.split('-')))
        return country, date(year, month, day)
    except ValueError:
        return country, None
    except:
        raise


def date2float(dt):
    origin = date(dt.year, 1, 1)
    td = (dt-origin).days
    return dt.year + td/365.25

def main():
    # load country to region mapping
    country2region = {}
    with open('countries.csv') as handle:
        for line in handle:
            try:
                country, region = line.strip().split(',')
            except:
                print(line)
                raise
            country2region.update({country: region})

    # tn93 -o data/gisaid.tn93.csv data/gisaid-filtered.fa
    handle = open('data/gisaid.tn93.csv')
    _ = next(handle)

    print("Filter sequences that are too distant")
    min_dists = {}
    unknown_country = []
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
        for country in unknown_country:
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


    print("building graph from nodes to find clusters")
    handle.seek(0)  # reset TN93 file
    _ = next(handle)
    G = nx.Graph()
    for line in handle:
        id1, id2, dist = line.strip().split(',')
        dist = float(dist)
        #if id1 in discard or id2 in discard:
        #    continue
        for node in [id1, id2]:
            if node not in G:
                country, coldate = parse_label(node)
                G.add_node(node, country=country, coldate=coldate)
        if dist < 1e-09:
            # some very small distances reported - mixtures?
            G.add_edge(id1, id2)

    components = list(nx.connected_components(G))
    print("graph comprises {} sequences and "
          "{} components".format(len(G), len(components)))

    # generate cluster information
    infofile = open('data/clusters.info.csv', 'w')
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
        clusters.update({label: coldate})

        # write cluster contents to info file
        for coldate, country, node in intermed:
            region = country2region[country]
            infofile.write('{},{},{},{},{}\n'.format(
                label, node, coldate, region, country
            ))


    outfile = open('data/clusters.fa', 'w')
    datefile = open('data/clusters.dates.csv', 'w')
    datefile.write('name,date\n')
    for h, s in iter_fasta(open('data/gisaid-filtered.fa')):
        if h in clusters:
            accession = h.split('|')[1]
            outfile.write(">{}\n{}\n".format(accession, s.replace('?', 'N')))
            datefile.write('{},{}\n'.format(
                accession, clusters[h]) #date2float(clusters[h]))
            )
    outfile.close()

# pass outputs to fasttree2 and treetime
# fasttree2 -nt < clusters.fa > clusters.ft2.nwk
# python3 prune-long-tips.py
# treetime --tree data/clusters.pruned.nwk --aln data/clusters.fa --dates data/clusters.dates.csv
# python3 parse-nexus.py
if __name__ == "__main__":
    main()

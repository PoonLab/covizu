import networkx as nx
from datetime import date
from csv import DictWriter
import sys

G = nx.Graph()

handle = open('data/gisaid-aligned.nopango.tn93.csv')
_ = next(handle)

# generate equivalence graph
countries = []
for line in handle:
    id1, id2, dist = line.strip().split(',')
    co1 = id1.split('/')[1]
    if co1 not in countries:
        countries.append(co1)
    co2 = id2.split('/')[1]
    if co2 not in countries:
        countries.append(co2)

    if id1 not in G:
        G.add_node(id1)
    if id2 not in G:
        G.add_node(id2)

    if float(dist) == 0:
        G.add_edge(id1, id2)

#clusters = [x for x in nx.connected_components(G) if len(x) > 1]

clust_file = DictWriter(open('clusters.csv', 'w'),
                        fieldnames=[
                            'label', 'size', 'mindate', 'maxdate',
                            'meandate'] + countries
                        )
clust_file.writeheader()

node2cluster = {}
for cluster in nx.connected_components(G):
    label = cluster.pop()  # arbitrary member
    node2cluster.update({label: label})

    # prepare containers
    country_counts = dict([(country, 0) for country in countries])
    dates = []

    for node in cluster:
        node2cluster.update({node: label})
        info, epi_id, ymd = node.split('|')

        country = info.split('/')[1]
        country_counts[country] += 1

        try:
            year, month, day = list(map(int, ymd.split('-')))
            dates.append(date(year, month, day).toordinal())
        except:
            # do not include missing date
            pass

    mindate = date.fromordinal(min(dates)) if dates else 'NA'
    maxdate = date.fromordinal(max(dates)) if dates else 'NA'
    meandate = date.fromordinal(round(sum(dates)/len(dates))) if dates else 'NA'
    # write cluster info to file
    row = dict(
        label=label, size=len(cluster),
        mindate=mindate, maxdate=maxdate, meandate=meandate
    )
    row.update(country_counts)
    clust_file.writerow(row)


sys.exit()

# generate graph collapsing identical sequences
handle.seek(0)  # reset file stream
_ = next(handle)
for line in handle:
    id1, id2, dist = line.strip().split(',')
    dist = float(dist)
    if dist == 0:
        continue
    c1 = node2cluster[id1]
    c2 = node2cluster[id2]

import networkx as nx
from networkx.algorithms import tree
from datetime import date
from csv import DictWriter
import math
import sys

G = nx.Graph()

# FIXME: we should exclude sequences with too many ambiguous nts
# tn93 -o gisaid-aligned.human.tn93.csv gisaid-aligned.human.fa
handle = open('data/gisaid-aligned.human.tn93.csv')
_ = next(handle)

# generate equivalence graph
print("building graph from nodes to find clusters")
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

    if float(dist) < 1e-6:
        # some very small distances reported - mixtures?
        G.add_edge(id1, id2)

components = list(nx.connected_components(G))
#sys.exit()

print("graph comprises {} sequences".format(len(G)))
#clusters = [x for x in nx.connected_components(G) if len(x) > 1]

clust_file = DictWriter(open('clusters.csv', 'w'),
                        fieldnames=[
                            'label', 'size', 'mindate', 'maxdate',
                            'meandate'] + countries
                        )
clust_file.writeheader()

def parse_label(label):
    info, epi_id, ymd = label.split('|')
    country = info.split('/')[1]
    try:
        year, month, day = list(map(int, ymd.split('-')))
        return country, date(year, month, day).toordinal()
    except ValueError:
        return country, None
    except:
        raise

print("annotating clusters")
node2cluster = {}
clusters = {}
canada = {}

for cluster in components:
    label = cluster.pop()  # arbitrary member
    clusters.update({label: len(cluster)+1})
    node2cluster.update({label: label})

    # prepare containers
    country_counts = dict([(country, 0) for country in countries])
    dates = []

    # parse label
    country, coldate = parse_label(label)
    country_counts[country] += 1
    if coldate:
        dates.append(coldate)

    for node in cluster:
        node2cluster.update({node: label})
        country, coldate = parse_label(node)

        country_counts[country] += 1
        if coldate:
            dates.append(coldate)

    mindate = date.fromordinal(min(dates)) if dates else 'NA'
    maxdate = date.fromordinal(max(dates)) if dates else 'NA'
    meandate = date.fromordinal(round(sum(dates)/len(dates))) if \
        dates else 'NA'

    # write cluster info to file
    row = dict(
        label=label, size=len(cluster)+1,
        mindate=mindate, maxdate=maxdate, meandate=meandate
    )
    row.update(country_counts)
    clust_file.writerow(row)

    # cache Canada status
    canada.update({label: country_counts['Canada'] /
                          (len(cluster)+1.)})


# generate graph collapsing identical sequences
# FIXME: this number is stochastic
print("building graph from {} clusters".format(len(clusters)))

G2 = nx.Graph()
handle.seek(0)  # reset file stream
_ = next(handle)
for line in handle:
    id1, id2, dist = line.strip().split(',')
    if id1 not in clusters or id2 not in clusters:
        # look exclusively at representative sequences
        continue

    if id1 not in G2:
        G2.add_node(id1, size=clusters[id1])
    if id2 not in G2:
        G2.add_node(id2, size=clusters[id2])

    if float(dist) < 1e-6:
        # this should never happen - problem with nx.connected_components?
        #G2.remove_node(id2)
        #nx.set_node_attributes(G2, {id1: clusters[id1]+clusters[id2]}, 'size')
        print("NOOO!")
        sys.exit()

    G2.add_edge(id1, id2, weight=float(dist))


# generate minimum spanning tree
print("generating minimum spanning tree")
mst = tree.minimum_spanning_tree(G2, weight='weight')

# write to dot file
print("generating DOT file")
dotfile = open('mst.dot', 'w')
dotfile.write('graph {\n')
dotfile.write('\tnode [shape="circle" style="filled" label=""];\n')

for node, size in mst.nodes(data='size'):
    dotfile.write('\t"{}" [width={} fillcolor="{}"];\n'.format(
        node,
        0.05 * math.sqrt(size)+0.05,
        #'red' if canada[node] else 'white'
        '#ff%02x%02x' % (round(255*(1-canada[node])),
                         round(255*(1-canada[node])))
    ))

# mean distance is 0.00035
for id1, id2, dist in mst.edges(data='weight'):
    dotfile.write('\t"{}" -- "{}" [len={}];\n'.format(
        id1, id2, 0.4)#dist/0.0003)
    )

dotfile.write("}\n")
dotfile.close()
# dot -Kneato -Tpdf mst.dot > mst.pdf

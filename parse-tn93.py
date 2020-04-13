import networkx as nx
from networkx.algorithms import tree
from datetime import date
from csv import DictWriter
import math
import sys


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


G = nx.Graph()

# FIXME: we should exclude sequences with too many ambiguous nts
# tn93 -o gisaid-aligned.human.tn93.csv gisaid-aligned.human.fa
handle = open('data/gisaid-aligned.human.tn93.csv')
_ = next(handle)

# generate equivalence graph
print("building graph from nodes to find clusters")
#countries = []
for line in handle:
    id1, id2, dist = line.strip().split(',')

    co1 = id1.split('/')[1]
    co2 = id2.split('/')[1]

    if id1 not in G:
        G.add_node(id1)
    if id2 not in G:
        G.add_node(id2)

    if float(dist) < 1e-6:
        # some very small distances reported - mixtures?
        G.add_edge(id1, id2)

components = list(nx.connected_components(G))


print("graph comprises {} sequences and "
      "{} components".format(len(G), len(components)))
#clusters = [x for x in nx.connected_components(G) if len(x) > 1]


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


# generate cluster information
clusters = {}
for cluster in components:
    intermed = [(parse_label(node)[1], node) for node in cluster]
    intermed = [(x, y) for x, y in intermed if x is not None]
    intermed.sort()  # increasing order of collection dates

    if len(intermed) == 0:
        #print('Omitting cluster of {} with no coldates'.format(
        #    len(cluster)
        #))
        continue

    label = intermed[0][1]  # label cluster by earliest case
    clusters.update({label: {}})
    for node in cluster:
        country, coldate = parse_label(node)
        region = country2region[country]

        if coldate is None:
            #print("Skipping node {} with no coldate".format(node))
            continue

        if coldate not in clusters[label]:
            clusters[label].update({coldate: []})
        clusters[label][coldate].append(region)


# generate graph collapsing identical sequences
print("building cluster graph".format(len(clusters)))

G2 = nx.Graph()
handle.seek(0)  # reset file stream
_ = next(handle)
for line in handle:
    id1, id2, dist = line.strip().split(',')
    if id1 not in clusters or id2 not in clusters:
        # look exclusively at representative sequences
        continue

    if id1 not in G2:
        G2.add_node(id1)
    if id2 not in G2:
        G2.add_node(id2)

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
edgelist = {}
for n1, n2, dist in mst.edges(data='weight'):
    if n1 not in edgelist:
        edgelist.update({n1: {}})
    edgelist[n1].update({n2: dist})

    if n2 not in edgelist:
        edgelist.update({n2: {}})
    edgelist[n2].update({n1: dist})


# traverse MST from earliest label (cluster)
# search edge list for children of root node and recurse
def traversal(node, parent, edgelist, history=[]):
    history.append(node)
    yield (node, parent)

    children = [child for child, _ in edgelist[node].items()
                if child not in history]

    # sort child nodes by earliest collection dates
    intermed = [(min(list(clusters[child].keys())), child)
                for child in children]
    intermed.sort()

    for _, child in intermed:
        for obj in traversal(child, node, edgelist, history=history):
            yield obj


# root the MST on the cluster with the most collection dates
intermed = [(len(clusters[node]), node) for node in clusters.keys()]
intermed.sort(reverse=True)
root = intermed[0][1]


# write clusters out to file
#dotfile = open('clusters.dot', 'w')
#dotfile.write('digraph {\n')
outfile = open('clusters.csv', 'w')

for child, parent in traversal(root, None, edgelist):
    #if parent is None:
    #    continue
    #dotfile.write('\t"{}"->"{}";\n'.format(
    #    parent.split('|')[0],
    #    child.split('|')[0]
    #))
    cases = list(clusters[child].items())
    cases.sort()
    for coldate, counts in cases:
        outfile.write(','.join([
            child,
            str(parent),
            str(coldate),
            '"{}"'.format(','.join(counts))
        ]))
        outfile.write('\n')

outfile.close()
#dotfile.write('}\n')
#dotfile.close()

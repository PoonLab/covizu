"""
Scratch file of minimum spanning tree code
"""

import networkx as nx

# load cluster data
nodes = {}
with open('data/clusters.info.csv') as f:
    _ = next(f)
    for line in f:
        label, node, coldate, region, country = line.strip().split(',')
        if label not in nodes:
            nodes.update({label: {
                'coldates': [],
                'regions': dict([(r, 0) for r in regions])
            }})
        nodes[label]['coldates'].append(parse_isodate(coldate))

# populate graph with nodes
G = nx.Graph()
for node, ndata in nodes.items():
    min_date =

# load cluster TN93
edges = []
with open('data/clusters.tn93.csv') as f:
    _ = next(f)
    for line in f:
        id1, id2, dist = line.strip().split(',')
        d = float(dist)
        if d < cutoff:
            edges.append((id1, id2, d))


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
#intermed = [(len(clusters[node]), node) for node in clusters.keys()]
#intermed.sort(reverse=True)
#root = intermed[0][1]

# root the MST on the cluster with the earliest collection date
intermed = [(min(list(cases.keys())), node) for node, cases in clusters.items()]
intermed.sort()
root = intermed[0][1]

# write clusters out to file
#dotfile = open('clusters.dot', 'w')
#dotfile.write('digraph {\n')
outfile = open('clusters.csv', 'w')
outfile.write("label,parent,coldate,country\n")

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

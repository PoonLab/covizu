import networkx as nx
from datetime import date

regions = ['Canada', 'USA', 'China', 'Asia', 'SAmerica', 'Europe',
           'Australia', 'Africa']

cutoff = 2e-4

def parse_isodate(iso):
    year, month, day = tuple(map(int, iso.split('-')))
    return date(year, month, day)

G = nx.Graph()

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

# load cluster TN93
edges = []
with open('data/clusters.tn93.csv') as f:
    _ = next(f)
    for line in f:
        id1, id2, dist = line.strip().split(',')
        d = float(dist)
        if d < cutoff:
            edges.append((id1, id2, d))

with open('data/clusters.tn93.dot', 'w') as f:
    f.write('graph {\n')
    f.write('  node [label="" shape="circle" style="filled"];\n')
    for node, nd in nodes.items():
        f.write('  "{}";\n'.format(node))

    for id1, id2, dist in edges:
        f.write('  "{}"--"{}";\n'.format(id1, id2))

    f.write('}\n')

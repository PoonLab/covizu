from datetime import date

def parse_label(label):
    info, epi_id, ymd = label.split('|')
    country = info.split('/')[1]
    ymd = list(map(int, ymd.split('-')))
    if len(ymd) == 3:
        year, month, day = ymd
        return country, date(year, month, day)
    elif len(ymd) == 2:
        year, month = ymd
        return country, date(year, month, 15)
    else:
        return country, None

country2region = {}
with open('countries.csv') as handle:
    for line in handle:
        try:
            country, region = line.strip().split(',')
        except:
            print(line)
            raise
        country2region.update({country: region})

cutoff = 2e-4

nodes = {}
edges = []

with open('data/gisaid.tn93.csv') as f:
    _ = next(f)
    for line in f:
        id1, id2, dist = line.strip().split(',')
        for node in [id1, id2]:
            if node not in nodes:
                country, coldate = parse_label(node)
                nodes.update({node: {
                    'coldate': coldate,
                    'region': country2region[country],
                    'country': country
                }})
        dist = float(dist)
        if dist < cutoff:
            edges.append((id1, id2, dist))

with open('data/gisaid.tn93.dot', 'w') as f:
    f.write('graph {\n')
    f.write('  node [label="" shape="circle" style="filled"];\n')
    for node, nd in nodes.items():
        f.write('  "{}";\n'.format(node))

    for id1, id2, dist in edges:
        f.write('  "{}"--"{}";\n'.format(id1, id2))

    f.write('};\n')

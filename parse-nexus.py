import re
import argparse
from clustering import date2float
from datetime import date
import sys
from Bio import Phylo
from io import StringIO

DATE_TOL = 0.1

parser = argparse.ArgumentParser(
    description = "Use regular expressions to extract comment fields "
                  "from NEXUS output of TreeTime and write to a "
                  "separate CSV file.  Remove problematic tips."
)
parser.add_argument('infile', type=argparse.FileType('r'),
                    help="input, TreeTime NEXUS output file")
parser.add_argument('csvfile', type=argparse.FileType('w'),
                    help="output, CSV file with node date estimates")
parser.add_argument('outfile', type=argparse.FileType('w'),
                    help="output, cleaned Newick file")
args = parser.parse_args()

handle = open('data/clusters.info.csv')
coldates = {}
for line in handle:
    _, node, dt, _, _ = line.strip().split(',')
    year, month, day = tuple(map(int, dt.split('-')))
    coldate = date(year, month, day)
    accession = node.split('|')[1]
    coldates.update({accession: date2float(coldate)})

# extract comment fields and store date estimates
pat = re.compile('([^)(,:]+):([0-9]+\.[0-9]+)\[[^d]+date=([0-9]+\.[0-9]+)\]')

# extract date estimates and internal node names
args.csvfile.write('node_name,branch_length,date_est\n')
remove = []
for line in args.infile:
    for m in pat.finditer(line):
        node_name, branch_length, date_est = m.groups()
        coldate = coldates.get(node_name, None)
        if coldate and (float(date_est) - coldate) > DATE_TOL:
            sys.stdout.write('removing {}:  {:0.3f} < {}\n'.format(
                node_name, coldate, date_est
            ))
            sys.stdout.flush()
            remove.append(node_name)

        args.csvfile.write('{},{},{}\n'.format(*m.groups()))
args.csvfile.close()

# second pass to excise all comment fields
pat = re.compile('\[&U\]|\[&mutations="[^"]*",date=[0-9]+\.[0-9]+\]')
args.infile.seek(0)
nexus = ''
for line in args.infile:
    nexus += pat.sub('', line)

# read in tree to prune problematic tips
phy = Phylo.read(StringIO(nexus), format='nexus')

for node_name in remove:
    phy.prune(node_name)

for node in phy.get_terminals():
    node.comment = None

for node in phy.get_nonterminals():
    node.name = node.confidence
    node.confidence = None
    node.comment = None

Phylo.write(phy, file=args.outfile, format='newick')

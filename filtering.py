from gotoh2 import iter_fasta
import re

MAX_PROP_N = 0.05
MIN_SEQ_LENGTH = 29000

#pat = re.compile('pangolin|canine|bat')
pat = re.compile('^[^/]+/[a-z]')
outfile = open('data/gisaid-filtered.fa', 'w')

n_non_human = 0
n_problem = 0
n_short = 0

for h, s in iter_fasta(open('data/gisaid-aligned.fa')):
    if pat.findall(h):
        #print("omitting non-human entry {}".format(h))
        n_non_human += 1
        continue
    if s.count('?') / float(len(s)) > MAX_PROP_N:
        #print("omitting entry {} with too many ?s".format(h))
        n_problem += 1
        continue
    if len(s.replace('-', '')) < MIN_SEQ_LENGTH:
        n_short += 1
        continue

    outfile.write('>{}\n{}\n'.format(h, s))

outfile.close()
print("Discarded {} non-human sequences.".format(n_non_human))
print("Discarded {} problematic sequences with prop N > {}.".format(
      n_problem, MAX_PROP_N))
print("Discarded {} sequences of length < {}".format(
    n_short, MIN_SEQ_LENGTH
))

from gotoh2 import iter_fasta
import re

#pat = re.compile('pangolin|canine|bat')
pat = re.compile('^[^/]+/[a-z]')
outfile = open('gisaid-aligned.human.fa', 'w')

for h, s in iter_fasta(open('data/gisaid-aligned.fa')):
    if pat.findall(h):
        print("omitting non-human entry {}".format(h))
        continue
    if s.count('?') / float(len(s)) > 0.1:
        print("omitting entry {} with too many ?s".format(h))
        continue

    outfile.write('>{}\n{}\n'.format(h, s))

outfile.close()

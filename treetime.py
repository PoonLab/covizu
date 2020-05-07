import argparse
import sys
from csv import DictReader
from gotoh2 import iter_fasta


def parse_clusters(infile):
    """
    :param infile: cluster info file from clustering.py
    :return: dict
    """
    clusters = {}
    reader = DictReader(open(infile))
    for row in reader:
        cluster = row['cluster']
        if cluster not in clusters:
            clusters.update({cluster: {
                'count': 0, 'coldate': row['coldate']
            }})
        clusters[cluster]['count'] += 1
    return clusters


def write_inputs(clusters, n_output, fasta_in, fasta_out, dates_out):
    """
    Export cluster sequences for TimeTree analysis
    :param n_output:  int, number of cluster sequences to output
    :param clusters:  dict, returned from write_info()
    :param fasta_in:  input, path to FASTA containing genome sequences
    :param fasta_out: output, path to write FASTA of cluster sequences
    :param dates_out: output, path to write cluster dates for TreeTime
    """
    # filter the largest clusters
    intermed = [(v['count'], k) for k, v in clusters.items()]
    intermed.sort(reverse=True)
    keys = [key for count, key in intermed[:n_output]]
    print(keys)

    # open file streams
    outfile = open(fasta_out, 'w')
    datefile = open(dates_out, 'w')
    datefile.write('name,date\n')

    for h, s in iter_fasta(open(fasta_in)):
        if h not in keys:
            continue
        accession = h.split('|')[1]
        outfile.write(">{}\n{}\n".format(accession, s.replace('?', 'N')))
        datefile.write('{},{}\n'.format(accession, clusters[h]['coldate']))

    outfile.close()

# pass outputs to fasttree2 and treetime
# fasttree2 -nt < clusters.fa > clusters.ft2.nwk
# python3 prune-long-tips.py
# treetime --tree data/clusters.pruned.nwk --aln data/clusters.fa --dates data/clusters.dates.csv
# python3 parse-nexus.py


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--info', default='data/clusters.info.csv',
                        help='Path to CSV file with cluster information.')
    parser.add_argument('--infasta', default='data/clusters.fa',
                        help='Path to FASTA file with aligned GISAID data.')
    parser.add_argument('--outfasta', default='data/treetime.fa',
                        help='Path to file to write FASTA output.')
    parser.add_argument('--dates', default='data/treetime.csv',
                        help='Path to file to write dates in CSV format '
                             '(for TreeTime analysis).')
    parser.add_argument('-n', default=20,
                        help='Number of sequences to export to TreeTime.')
    return parser.parse_args()


if __name__ == '__main__':
    def callback(msg):
        sys.stdout.write(msg+'\n')
        sys.stdout.flush()

    args = parse_args()
    clusters = parse_clusters(args.info)
    print(clusters)
    write_inputs(clusters, n_output=args.n, fasta_in=args.infasta,
                 fasta_out=args.outfasta, dates_out=args.dates)

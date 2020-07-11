import re
import argparse
import sqlite3


def get_aligned(database):
    """
    Connect to sqlite3 database and stream header, aligned genome tuples.
    TODO: allow user to limit query to range of sample collection dates?

    :param database:  str, path to sqlite3 database
    :yield:  header, sequence
    """
    conn = sqlite3.connect(database, check_same_thread=False)
    cur = conn.cursor()
    data = cur.execute('SELECT `header`, `aligned` FROM SEQUENCES;').fetchall()
    for h, s in data:
        yield h, s


def filter_gisaid(database, outfile, trim_left=0, trim_right=0,
                  max_prop_n=0.05, minlen=29000):
    """
    Filter FASTA file for partial and non-human SARS-COV-2 genome sequences.

    :param fasta_file:  open file stream to GISAID FASTA file
    :param outfile:  open file stream to write filtered FASTA file
    :param trim_left:  int, number of bases to drop from left
    :param trim_right:  int, number of bases to drop from right
    :param max_prop_n:  float, maximum proportion of N's (ambiguous bases)
                        tolerated per genome
    :param minlen:  int, minimum tolerated sequence length

    :return:  dict, containing lists of headers for rejected genomes
    """
    # lower-case label in place of country identifies non-human samples
    pat = re.compile('^[^/]+/[a-z]')
    pat2 = re.compile("^[HhCcOoVv]+-19/[A-Z][^/]+/[^/]+/[0-9-]+\|[^|]+\|[0-9]{4}-[0-9]+-[0-9]+")
    pat3 = re.compile('^-*')
    pat4 = re.compile('-*$')

    accessions = {}
    discards = {'nonhuman': [], 'ambiguous': [], 'short': [],
                'duplicates': [], 'mangled header': []}

    for h, s in get_aligned(database):
        if not type(h)==str or not type(s)==str:
            print("Error: entry {} not string type: sequence {}".format(h, s))
            continue

        if pat.findall(h):
            discards['nonhuman'].append(h)
            continue
        
        if len(s) < minlen:
            discards['short'].append(h)
            continue

        # apply sequence trims
        seq = s[trim_left:(-trim_right)]

        # this is conservative - all internal gaps are interpreted as deletions
        gap_prefix = len(pat3.findall(seq)[0])
        gap_suffix = len(pat4.findall(seq)[0])
        seqlen = len(seq) - gap_prefix - gap_suffix

        n_ambig = seq.count('?') + seq.count('N') + gap_prefix + gap_suffix
        if n_ambig / float(len(seq)) > max_prop_n:
            discards['ambiguous'].append(h)
            continue

        if pat2.search(h) is None:
            discards['mangled header'].append(h)
            continue

        desc, accn, coldate = h.split('|')
        if accn in accessions:
            discards['duplicates'].append(h)
            continue
        accessions.update({accn: desc})

        # write genome to output file
        _ = outfile.write('>{}\n{}\n'.format(h, seq))

    return discards


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter aligned sequences for problematic entries."
    )

    parser.add_argument('-i', '--database', type=str, default='data/gsaid.db',
                        help='input, path to sqlite3 database.')

    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=open('data/gisaid-filtered.fa', 'w'),
                        help='output, path to write filtered FASTA file')

    parser.add_argument('-p', '--maxpropN', type=float, default=0.05,
                        help='option, maximum tolerance for proportion of Ns '
                             '(ambiguous base calls) - defaults to 0.05 (5%)')

    parser.add_argument('-L', '--minlen', type=int, default=29000,
                        help='option, minimum genome sequence length, default 29000')

    parser.add_argument('--trim_left', type=int, default=53,
                        help='option, remove N bases from the left')

    parser.add_argument('--trim_right', type=int, default=93,
                        help='option, remove N bases from the right')

    parser.add_argument('--verbose', action='store_true',
                        help='option, print discarded genome headers')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    discards = filter_gisaid(
        database=args.database, outfile=args.outfile, trim_left=args.trim_left,
        trim_right=args.trim_right, max_prop_n=args.maxpropN, minlen=args.minlen
    )

    print("Discarded {} non-human sequences.".format(len(discards['nonhuman'])))
    if args.verbose:
        for h in discards['nonhuman']:
            print('  {}'.format(h))

    print("Discarded {} problematic sequences with prop N > {}.".format(
        len(discards['ambiguous']), args.maxpropN
    ))
    if args.verbose:
        for h in discards['ambiguous']:
            print('  {}'.format(h))

    print("Discarded {} sequences of length < {}".format(
        len(discards['short']), args.minlen
    ))
    if args.verbose:
        for h in discards['short']:
            print('  {}'.format(h))

    print("Discarded {} sequences with duplicate accession numbers".format(
        len(discards['duplicates'])
    ))
    if args.verbose:
        for h in discards['duplicates']:
            print('  {}'.format(h))

    print("Discarded {} sequences with mangled headers / ambiguous sample dates".format(
        len(discards['mangled header'])
    ))
    if args.verbose:
        for h in discards['mangled header']:
            print('  {}'.format(h))

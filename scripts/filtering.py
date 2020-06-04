from gotoh2 import iter_fasta
import re
import argparse



def filter_gisaid(fasta_file, outfile, max_prop_n=0.05, minlen=29000):
    """
    Filter FASTA file for partial and non-human SARS-COV-2 genome sequences.

    :param fasta_file:  open file stream to GISAID FASTA file
    :param outfile:  open file stream to write filtered FASTA file
    :param max_prop_n:  maximum proportion of N's (ambiguous bases) tolerated per genome
    :param minlen:  minimum tolerated sequence length

    :return:  dict, containing lists of headers for rejected genomes
    """
    # lower-case label in place of country identifies non-human samples
    pat = re.compile('^[^/]+/[a-z]')
    pat2 = re.compile("^[HhCcOoVv]+-19/[A-Z][^/]+/[^/]+/[0-9-]+\|[^|]+\|[0-9]{4}-[0-9]+-[0-9]+")

    accessions = {}
    discards = {'nonhuman': [], 'ambiguous': [], 'short': [],
                'duplicates': [], 'mangled header': []}

    for h, s in iter_fasta(fasta_file):
        if pat.findall(h):
            discards['nonhuman'].append(h)
            continue
        if s.count('?') / float(len(s)) > max_prop_n:
            discards['ambiguous'].append(h)
            continue
        if len(s.replace('-', '')) < minlen:
            discards['short'].append(h)
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
        outfile.write('>{}\n{}\n'.format(h, s))

    return discards


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter GISAID FASTA file for incomplete and non-human genome sequences."
    )
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=open('data/gisaid-aligned.fa'),
                        help='input, path to FASTA file with aligned GISAID genomes')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=open('data/gisaid-filtered.fa', 'w'),
                        help='output, path to write filtered FASTA file')
    parser.add_argument('-p', '--maxpropN', type=float, default=0.05,
                        help='option, maximum tolerance for proportion of Ns '
                             '(ambiguous base calls) - defaults to 0.05 (5%)')
    parser.add_argument('-L', '--minlen', type=int, default=29000,
                        help='option, minimum genome sequence length, default 29000')
    parser.add_argument('--verbose', action='store_true',
                        help='option, print discarded genome headers')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    discards = filter_gisaid(args.infile, args.outfile, args.maxpropN, args.minlen)

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
    for h in discards['duplicates']:
        print('  {}'.format(h))

    print("Discarded {} sequences with mangled headers".format(
        len(discards['mangled header'])
    ))
    for h in discards['mangled header']:
        print('  {}'.format(h))

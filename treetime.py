import argparse
from subprocess import Popen, PIPE, check_call
import json
from gotoh2 import iter_fasta
from tempfile import NamedTemporaryFile
import os
import math


def filter_fasta(fasta_file, json_file, cutoff=10):
    """
    Filter the variants FASTA file for genomes representing clusters
    identified by hierarchical clustering (hclust.R).  Add variants
    that are observed N > [cutoff] times in the data.

    :param fasta_file:  path to FASTA file containing cluster sequences
    :param json_file:  path to JSON file with cluster information
    :return:  dict, filtered header-sequence pairs
    """
    result = {}
    fasta = dict(list(iter_fasta(fasta_file)))
    clusters = json.load(json_file)
    for cluster in clusters:
        # record variant in cluster that is closest to root
        if type(cluster['nodes']) is list:
            # omit problematic cluster of one
            print(cluster['nodes'])
            continue

        header = list(cluster['nodes'].keys())[0]
        result.update({header: fasta[header]})

        # extract variants in cluster that have high counts
        major = [label for label, samples in
                 cluster['nodes'].items() if
                 len(samples) > cutoff and label != header]
        for label in major:
            result.update({label: fasta[label]})

    return result


def fasttree(fasta):
    """
    Wrapper for FastTree2, passing FASTA as stdin and capturing the
    resulting Newick tree string as stdout.
    :param fasta: dict, header: sequence pairs
    :return: str, Newick tree string
    """
    in_str = ''
    for h, s in fasta.items():
        accn = h.split('|')[1]
        in_str += '>{}\n{}\n'.format(accn, s)
    p = Popen(['fasttree2', '-nt', '-quote'], stdin=PIPE, stdout=PIPE)
    # TODO: exception handling with stderr?
    stdout, stderr = p.communicate(input=in_str.encode('utf-8'))
    return stdout.decode('utf-8')


def treetime(nwk, fasta, outdir, clock=None):
    """
    :param nwk: str, Newick tree string from fasttree()
    :param fasta: dict, header-sequence pairs
    :param outdir:  path to write output files
    :param clock: float, clock rate to constrain analysis - defaults
                  to None (no constraint)
    :return:  path to NEXUS output file
    """
    # extract dates from sequence headers
    datefile = NamedTemporaryFile('w', delete=False)
    datefile.write('name,date\n')
    alnfile = NamedTemporaryFile('w', delete=False)
    for h, s in fasta.items():
        # TreeTime seems to have trouble handling labels with spaces
        _, accn, coldate = h.split('|')
        datefile.write('{},{}\n'.format(accn, coldate))
        alnfile.write('>{}\n{}\n'.format(accn, s))
    datefile.close()
    alnfile.close()

    with NamedTemporaryFile('w', delete=False) as nwkfile:
        nwkfile.write(nwk.replace(' ', ''))

    call = ['treetime', '--tree', nwkfile.name,
            '--aln', alnfile.name, '--dates', datefile.name,
            '--outdir', outdir]
    if clock:
        call.extend(['--clock-rate', str(clock)])
    check_call(call)

    nexus_file = os.path.join(outdir, 'timetree.nexus')
    if not os.path.exists(nexus_file):
        print("Error: missing expected NEXUS output file {}".format(nexus_file))
        return None
    return nexus_file


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate inputs for TreeTime analysis."
    )
    parser.add_argument('--json', type=argparse.FileType('r'),
                        default=open('data/clusters.json'),
                        help='input, JSON file generated by hclust.R '
                             'identifying representative cluster '
                             'sequences')
    parser.add_argument('--fasta', type=argparse.FileType('r'),
                        default=open('data/variants.fa'),
                        help='input, FASTA file with unique variant '
                             'sequences')
    parser.add_argument('--mincount', type=int, default=math.inf,
                        help='optional, minimum count of variant to be '
                             'added to tree (default: infinity)')
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='optional, specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')
    parser.add_argument('--outdir', default='treetime/',
                        help='directory to write TreeTime output files')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    fasta = filter_fasta(args.fasta, args.json, cutoff=args.mincount)
    nwk = fasttree(fasta)
    treetime(nwk, fasta, outdir=args.outdir, clock=args.clock)

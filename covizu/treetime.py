import argparse
from subprocess import Popen, PIPE, check_call
import json
from tempfile import NamedTemporaryFile
import os
from io import StringIO

from Bio import Phylo

import covizu
from covizu.utils.seq_utils import *
from covizu.utils.db_utils import dump_raw_by_lineage, retrieve_seqs
from covizu.utils.progress_utils import Callback
from covizu.minimap2 import minimap2, encode_diffs


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
    fasta = dict([(h.split('|')[1], {'sequence': s, 'label': h}) for
                  h, s in iter_fasta(fasta_file)])
    clusters = json.load(json_file)
    for cluster in clusters:
        # record variant in cluster that is closest to root
        if type(cluster['nodes']) is list:
            # omit problematic cluster of one
            print(cluster['nodes'])
            continue

        # first entry is representative variant
        accn = list(cluster['nodes'].keys())[0]
        result.update({fasta[accn]['label']: fasta[accn]['sequence']})

        # extract other variants in cluster that have high counts
        major = [label for label, samples in
                 cluster['nodes'].items() if
                 len(samples) > cutoff and label != accn]
        for label in major:
            result.update({fasta[label]['label']: fasta[label]['sequence']})

    return result


def fasttree(fasta, binpath='fasttree2', seed=1):
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
    p = Popen([binpath, '-nt', '-quote', '-seed', str(seed)],
              stdin=PIPE, stdout=PIPE)
    # TODO: exception handling with stderr?
    stdout, stderr = p.communicate(input=in_str.encode('utf-8'))
    return stdout.decode('utf-8')


def treetime(nwk, fasta, outdir, binpath='treetime', clock=None, verbosity=1):
    """
    :param nwk: str, Newick tree string from fasttree()
    :param fasta: dict, header-sequence pairs
    :param outdir:  path to write output files
    :param clock: float, clock rate to constrain analysis - defaults
                  to None (no constraint)
    :param verbosity:  verbose level, defaults to 1
    :return:  path to NEXUS output file
    """
    # extract dates from sequence headers
    datefile = NamedTemporaryFile('w', prefix="cvz_tt_", delete=False)
    datefile.write('name,date\n')
    alnfile = NamedTemporaryFile('w', prefix="cvz_tt_", delete=False)
    for h, s in fasta.items():
        # TreeTime seems to have trouble handling labels with spaces
        _, accn, coldate = h.split('|')
        datefile.write('{},{}\n'.format(accn, coldate))
        alnfile.write('>{}\n{}\n'.format(accn, s))
    datefile.close()
    alnfile.close()

    with NamedTemporaryFile('w', prefix="cvz_tt_", delete=False) as nwkfile:
        nwkfile.write(nwk.replace(' ', ''))

    call = [binpath, '--tree', nwkfile.name,
            '--aln', alnfile.name, '--dates', datefile.name,
            '--outdir', outdir, '--verbose', str(verbosity)]
    if clock:
        call.extend(['--clock-rate', str(clock)])
    check_call(call)

    # clean up temporary files
    os.remove(datefile.name)
    os.remove(alnfile.name)
    os.remove(nwkfile.name)

    # return path to NEXUS file
    nexus_file = os.path.join(outdir, 'timetree.nexus')
    if not os.path.exists(nexus_file):
        print("Error: missing expected NEXUS output file {}".format(nexus_file))
        return None
    return nexus_file


def date2float(isodate):
    """ Convert ISO date string to float """
    year, month, day = map(int, isodate.split('-'))
    dt = date(year, month, day)
    origin = date(dt.year, 1, 1)
    td = (dt-origin).days
    return dt.year + td/365.25


def parse_nexus(nexus_file, fasta, date_tol):
    """
    Converting Treetime NEXUS output into Newick

    @param nexus_file:  str, path to write Newick tree string
    @param fasta:  dict, {header: seq} from filter_fasta()
    @param date_tol:  float, tolerance in tip date discordance
    """
    coldates = {}
    for h, _ in fasta.items():
        _, accn, coldate = h.split('|')
        coldates.update({accn: date2float(coldate)})

    # extract comment fields and store date estimates
    pat = re.compile('([^)(,:]+):([0-9]+\.[0-9]+)\[[^d]+date=([0-9]+\.[0-9]+)\]')

    # extract date estimates and internal node names
    remove = []
    with open(nexus_file) as handle:
        for line in handle:
            for m in pat.finditer(line):
                node_name, branch_length, date_est = m.groups()
                coldate = coldates.get(node_name, None)
                if coldate and abs(float(date_est) - coldate) > date_tol:
                    sys.stdout.write('removing {}:  {:0.3f} < {}\n'.format(
                        node_name, coldate, date_est
                    ))
                    sys.stdout.flush()
                    remove.append(node_name)

    # second pass to excise all comment fields
    pat = re.compile('\[&U\]|\[&mutations="[^"]*",date=[0-9]+\.[0-9]+\]')
    nexus = ''
    for line in open(nexus_file):
        nexus += pat.sub('', line)

    # read in tree to prune problematic tips
    phy = Phylo.read(StringIO(nexus), format='nexus')
    for node_name in remove:
        phy.prune(node_name)

    for node in phy.get_terminals():
        node.comment = None

    for node in phy.get_nonterminals():
        if node.name is None and node.confidence:
            node.name = node.confidence
            node.confidence = None
        node.comment = None

    Phylo.write(phy, file=nexus_file.replace('.nexus', '.nwk'),
                format='newick')


def retrieve_genomes(db="data/gsaid.db", stream=False, nthread=1, ref_file='data/MT291829.fa', misstol=300):
    """
    Query database for Pangolin lineages and then retrieve the earliest
    sampled genome sequence for each.  Export as FASTA for TreeTime analysis.
    :param db:  str, path to sqlite3 database
    :return:  list, (header, sequence) tuples
    """

    # load and parse reference genome
    with open(ref_file) as handle:
        _, refseq = convert_fasta(handle)[0]
    reflen = len(refseq)

    # allocate lists
    coldates = []
    lineages = []
    seqs = []

    # retrieve unaligned genomes from database
    for lineage, fasta_file in dump_raw_by_lineage(db):
        mm2 = minimap2(infile=fasta_file, nthread=nthread, stream=stream, ref=ref_file)
        gen = encode_diffs(mm2, reflen=reflen)
        for row in filter_outliers(gen):
            # exclude genomes too divergent from expectation
            if total_missing(row) > misstol:
                continue
            # take the earliest valid genome
            qname, _, _ = row
            _, coldate = parse_label(qname)
            break

        # update lists
        lineages.append(lineage)
        coldates.append(coldate)

        # reconstruct aligned sequence from feature vector
        seq = apply_features(row, refseq=refseq)
        seqs.append(seq)

        # clean up temporary files
        os.remove(fasta_file.name)

    # generate new headers in {name}|{accession}|{date} format expected by treetime()
    headers = map(lambda xy: '|{}|{}'.format(*xy), zip(lineages, coldates))
    return dict(zip(headers, seqs))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate inputs for TreeTime analysis."
    )
    parser.add_argument('--db', type=str, default='data/gsaid.db',
                        help='input, sqlite3 database')
    parser.add_argument('--ref', type=argparse.FileType('r'),
                        default=open(os.path.join(covizu.__path__[0]), "data/MT291829.fa"),
                        help="input, FASTA file with reference genome")
    parser.add_argument('--misstol', type=int, default=300,
                        help="optional, maximum tolerated number of missing bases per "
                             "genome (default 300).")
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='optional, specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')
    # FIXME: I think the next argument is deprecated
    parser.add_argument('--datetol', type=float, default=0.1,
                        help='optional, exclude tips from time-scaled tree '
                             'with high discordance between estimated and '
                             'known sample collection dates (year units,'
                             'default: 0.1)')
    parser.add_argument('--outdir', default='data/',
                        help='optional, directory to write TreeTime output files')
    parser.add_argument('--ft2bin', default='fasttree2',
                        help='optional, path to fasttree2 binary executable')
    parser.add_argument('--ttbin', default='treetime',
                        help='optional, path to treetime binary executable')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cb = Callback()

    cb.callback("Retrieving genomes")
    fasta = retrieve_genomes(args.db, ref_file=args.ref, misstol=args.misstol)

    cb.callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = fasttree(fasta, binpath=args.ft2bin)

    cb.callback("Reconstructing time-scaled tree with {}").format(args.ttbin)
    nexus_file = treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                          clock=args.clock)

    cb.callback("")
    parse_nexus(nexus_file, fasta, date_tol=args.datetol)

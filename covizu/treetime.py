import argparse
from subprocess import Popen, PIPE, check_call
import json
from tempfile import NamedTemporaryFile
import os
import sys
from io import StringIO
import re
import statistics  # Python 3.4+

from Bio import Phylo

import covizu
from covizu.utils.seq_utils import *
from covizu.utils.progress_utils import Callback


def fasttree(fasta, binpath='fasttree2', seed=1, gtr=True, collapse=True):
    """
    Wrapper for FastTree2, passing FASTA as stdin and capturing the
    resulting Newick tree string as stdout.

    :param fasta:  dict, header: sequence pairs
    :param binpath:  str, path to Fasttree2 binary executable
    :param seed:  int, initialize random number generator
    :param gtr:  if True, use generalized time reversible model; otherwise JC

    :return: str, Newick tree string
    """
    # convert dict to FASTA-formatted string
    in_str = ''
    for h, s in fasta.items():
        accn = h.split('|')[1]
        in_str += '>{}\n{}\n'.format(accn, s)

    cmd = [binpath, '-nt', '-quote', '-seed', str(seed)]
    if gtr:
        cmd.append('-gtr')

    p = Popen(cmd, stdin=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate(input=in_str.encode('utf-8'))
    nwk = stdout.decode('utf-8')

    # collapse low support nodes into polytomies
    if collapse:
        phy = Phylo.read(StringIO(nwk), format='newick')
        phy.collapse_all(lambda x: x.confidence is not None and
                                   x.confidence < 0.5)
        nwk = phy.format('newick')
    return nwk


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
            '--outdir', outdir, '--verbose', str(verbosity),
            '--plot-rtt', 'none',  # see issue #66
            '--clock-filter', '0',  # issue #245
            '--keep-polytomies'  # issue #339
    ]
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
    """ Convert ISO date string to float (years) """
    year, month, day = map(int, isodate.split('-'))
    dt = date(year, month, day)
    origin = date(dt.year, 1, 1)
    td = (dt-origin).days
    return dt.year + td/365.25


def parse_nexus(nexus_file, fasta, callback=None):
    """
    Converting Treetime NEXUS output into Newick

    @param nexus_file:  str, path to TreeTime NEXUS output
    @param fasta:  dict, {header: seq} from retrieve_genomes()
    @param callback:  function, optional callback
    """
    coldates = {}
    for h, _ in fasta.items():
        _, accn, coldate = h.split('|')
        coldates.update({accn: date2float(coldate)})

    # extract comment fields and store date estimates
    pat = re.compile('([^)(,:]+):(-*[0-9]+\.[0-9]+)\[[^d]+date=([0-9]+\.[0-9]+)\]')

    # extract date estimates and internal node names
    remove = []
    residuals = {}
    with open(nexus_file) as handle:
        for line in handle:
            for m in pat.finditer(line):
                node_name, branch_length, date_est = m.groups()
                coldate = coldates.get(node_name, None)  # internal nodes have no collection date
                if coldate:
                    # if estimated date is more recent (higher) than actual date,
                    # lineage has more substitutions than expected under molecular clock
                    residuals.update({node_name: float(date_est) - coldate})

    # second pass to excise all comment fields
    pat = re.compile('\[&U\]|\[&mutations="[^"]*",date=[0-9]+\.[0-9]+\]')
    nexus = ''
    for line in open(nexus_file):
        nexus += pat.sub('', line)

    # read in tree to prune problematic tips
    phy = Phylo.read(StringIO(nexus), format='nexus')
    for node_name in remove:
        phy.prune(node_name)

    # normalize residuals and append to tip labels
    rvals = residuals.values()
    rmean = statistics.mean(rvals)
    rstdev = statistics.stdev(rvals)
    for tip, resid in residuals.items():
        residuals[tip] = (resid-rmean) / rstdev

    for node in phy.get_terminals():
        node.comment = None

    for node in phy.get_nonterminals():
        if node.name is None and node.confidence:
            node.name = node.confidence
            node.confidence = None
        node.comment = None

    return phy, residuals


def retrieve_genomes(by_lineage, known_seqs, ref_file, earliest=True, callback=None):
    """
    Identify most recent sampled genome sequence for each Pangolin lineage.
    Export as FASTA for TreeTime analysis.

    :param by_lineage:  dict, return value from gisaid_utils::sort_by_lineage
    :param known_seqs:  dict, sequences used in Pango lineage designations
    :param ref_file:  str, path to FASTA file containing reference genome
    :param earliest:  bool, if False then use most recent genome as lineage representative
    :param callback:  optional, callback function

    :return:  dict, aligned genome keyed by header "|{lineage}|{coldate}"
    """
    # load and parse reference genome
    with open(ref_file) as handle:
        _, refseq = convert_fasta(handle)[0]

    # allocate lists
    coldates = []
    lineages = []
    seqs = []

    # retrieve unaligned genomes from database
    for lineage, records in by_lineage.items():
        # filter records for lineage-defining genomes
        curated = filter(
            lambda r: r['label'].replace(' ', '_')  # issue #313
                      in known_seqs,
            records)

        curated = list(curated)  # resolve generator
        if len(curated) == 0:
            if callback:
                callback("Error in retrieve_genomes(): no sequence names for lineage {} in designated "
                         "list; may need to update data/lineages.csv".format(lineage), level='WARN')
            curated = records

        intermed = [(r['coldate'], r['diffs'], r['missing']) for r in curated]
        intermed.sort(reverse=True)  # descending order
        coldate, diffs, missing = intermed[-1] if earliest else intermed[0]

        # update lists
        lineages.append(lineage)
        coldates.append(coldate)

        # reconstruct aligned sequence from feature vector
        seq = apply_features(diffs, missing=missing, refseq=refseq)
        seqs.append(seq)

    # generate new headers in {name}|{accession}|{date} format expected by treetime()
    headers = map(lambda xy: '|{}|{}'.format(*xy), zip(lineages, coldates))
    return dict(zip(headers, seqs))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate inputs for TreeTime analysis."
    )
    parser.add_argument('json', type=str, help='input, JSON produced by gisaid_utils.py')

    # issue #230
    path = covizu.__path__[0]
    if '.egg' in path:
        import pkg_resources
        ref_file = pkg_resources.resource_filename('covizu', 'data/NC_045512.fa')
    else:
        ref_file = os.path.join(path, "data", "NC_045512.fa")

    parser.add_argument('--ref', type=str, default=ref_file,
                        help="input, FASTA file with reference genome")
    parser.add_argument('--misstol', type=int, default=300,
                        help="optional, maximum tolerated number of missing bases per "
                             "genome (default 300).")
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='optional, specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')

    parser.add_argument('--outdir', default='data/',
                        help='optional, directory to write TreeTime output files')
    parser.add_argument('--ft2bin', default='fasttree2',
                        help='optional, path to fasttree2 binary executable')
    parser.add_argument('--ttbin', default='treetime',
                        help='optional, path to treetime binary executable')
    parser.add_argument('--lineages', type=str,
                        default=os.path.join(covizu.__path__[0], "data/pango-designation/lineages.csv"),
                        help="optional, path to CSV file containing Pango lineage designations.")

    parser.add_argument('--outfile', default='data/timetree.nwk',
                        help='output, path to write Newick tree string')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cb = Callback()

    cb.callback("Retrieving genomes")
    with open(args.json) as handle:
        by_lineage = json.load(handle)

    cb.callback("Parsing Pango lineage designations")
    handle = open(args.lineages)
    header = next(handle)
    if header != 'taxon,lineage\n':
        cb.callback("Error: {} does not contain expected header row 'taxon,lineage'".format(args.lineages))
        sys.exit()
    lineages = {}
    for line in handle:
        try:
            taxon, lineage = line.strip().split(',')
            if taxon and lineage:
                lineages.update({taxon: lineage})
            else:
                cb.callback("Warning '{}': taxon or lineage is missing".format(line), level='WARN')
        except:
            cb.callback("Warning: There is an issue with the line '{}' in lineages.csv".format(line), level='WARN')

    cb.callback("Identifying lineage representative genomes")
    fasta = retrieve_genomes(by_lineage, known_seqs=lineages, ref_file=args.ref, earliest=args.earliest,
                             callback=cb.callback)

    cb.callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = fasttree(fasta, binpath=args.ft2bin)

    cb.callback("Reconstructing time-scaled tree with {}".format(args.ttbin))
    nexus_file = treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                          clock=args.clock)

    timetree, residuals = parse_nexus(nexus_file, fasta)
    Phylo.write(timetree, file=args.outfile, format='newick')

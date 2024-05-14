"""what treetime does"""
import argparse
from subprocess import Popen, PIPE, check_call
import json
from tempfile import NamedTemporaryFile
import os
import sys
from io import StringIO
import re
import statistics  # Python 3.4+

from datetime import date
from Bio import Phylo
import pkg_resources
import covizu
from covizu.utils.seq_utils import convert_fasta, apply_features
from covizu.utils.progress_utils import Callback
# from covizu.utils.batch_utils import unpack_records
import covizu.utils.batch_utils


def fasttree(in_fasta, binpath='fasttree2', seed=1, gtr=True, collapse=True):
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
    for head, seq in in_fasta.items():
        accn = head.split('|')[1]
        in_str += f'>{accn}\n{seq}\n'

    cmd = [binpath, '-nt', '-quote', '-seed', str(seed)]
    if gtr:
        cmd.append('-gtr')

    with Popen(cmd, stdin=PIPE, stdout=PIPE) as p_file:
        stdout = p_file.communicate(input=in_str.encode('utf-8'))[0]
        newick = stdout.decode('utf-8')

        # parse Newick tree string from output
        phy = Phylo.read(StringIO(newick), format='newick')
        phy.root_with_outgroup('reference')  # issue #396
        phy.prune('reference')
        if collapse:
            # collapse low support nodes into polytomies
            phy.collapse_all(
                lambda x: x.confidence is not None and x.confidence < 0.5)
        newick = phy.format('newick')
        return newick


def treetime(newick, in_fasta, outdir, binpath='treetime', clock=None, verbosity=1):
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
    with NamedTemporaryFile('w', prefix="cvz_tt_", delete=False) as datefile:
        datefile.write('name,date\n')
        with NamedTemporaryFile('w', prefix="cvz_tt_", delete=False) as alnfile:
            for head, sequence in in_fasta.items():
                # TreeTime seems to have trouble handling labels with spaces
                _, accn, coldate = head.split('|')
                datefile.write(f'{accn},{coldate}\n')
                alnfile.write(f'>{accn}\n{sequence}\n')

    with NamedTemporaryFile('w', prefix="cvz_tt_", delete=False) as nwkfile:
        nwkfile.write(newick.replace(' ', ''))

    call = [binpath, '--tree', nwkfile.name,
            '--aln', alnfile.name, '--dates', datefile.name,
            '--outdir', outdir, '--verbose', str(verbosity),
            '--plot-rtt', 'none',  # see issue #66
            '--clock-filter', '0',  # issue #245
            '--keep-polytomies',  # issue #339
            '--keep-root'  # issue 396
            ]
    if clock:
        call.extend(['--clock-rate', str(clock)])
    check_call(call)

    # clean up temporary files
    os.remove(datefile.name)
    os.remove(alnfile.name)
    os.remove(nwkfile.name)

    # return path to NEXUS file
    nex_file = os.path.join(outdir, 'timetree.nexus')
    if not os.path.exists(nex_file):
        print(f"Error: missing expected NEXUS output file {nex_file}")
        return None
    return nex_file


def date2float(isodate):
    """ Convert ISO date string to float (years) """
    year, month, day = map(int, isodate.split('-'))
    date_time = date(year, month, day)
    origin = date(date_time.year, 1, 1)
    today = (date_time - origin).days
    return date_time.year + today / 365.25


def parse_nexus(in_nexus_file, in_fasta, callback=None):
    """
    Converting Treetime NEXUS output into Newick

    @param nexus_file:  str, path to TreeTime NEXUS output
    @param fasta:  dict, {header: seq} from retrieve_genomes()
    @param callback:  function, optional callback
    """
    coldates = {}
    for head, _ in in_fasta.items():
        if 'reference' in head:
            continue
        _, accn, coldate = head.split('|')
        coldates.update({accn: date2float(coldate)})

    # extract date estimates and internal node names
    result, to_rem = extract_nexus(in_nexus_file, coldates)

    phy = prune_nexus(in_nexus_file, to_rem)

    normalize_residule(result, callback)

    return edit_py_node(phy), result

def extract_nexus(in_nexus_file, coldates):
    """extract comment fields and store date estimates"""
    remove = []
    res = {}
    pat = re.compile(
        '([^)(,:]+):(-*[0-9]+\\.[0-9]+)\\[[^d]+date=([0-9]+\\.[0-9]+)\\]')
    with open(in_nexus_file, encoding='utf-8') as nex:
        for segment in nex:
            for matches in pat.finditer(segment):
                node_name, date_est = matches.groups()[0], matches.groups()[2]
                # internal nodes have no collection date
                coldate = coldates.get(node_name, None)
                if coldate:
                    # if estimated date is more recent (higher) than actual date,
                    # lineage has more substitutions than expected under
                    # molecular clock
                    res.update({node_name: float(date_est) - coldate})
    return res, remove

def prune_nexus(in_nex_file, remove):
    """second pass to excise all comment fields"""
    pat = re.compile('\\[&U\\]|\\[&mutations="[^"]*",date=[0-9]+\\.[0-9]+\\]')
    nex = ''

    with open(in_nex_file, encoding='utf-8') as nex:
        for segment in nex:
            nex += pat.sub('', segment)

    # read in tree to prune problematic tips
    phy = Phylo.read(StringIO(nex), format='nexus')
    for node_name in remove:
        phy.prune(node_name)

    return phy

def normalize_residule(in_res, in_call):
    """normalize residuals and append to tip labels"""
    rvals = in_res.values()
    try:
        rmean = statistics.mean(rvals)
        rstdev = statistics.stdev(rvals)
    except statistics.StatisticsError:
        in_call("Provided records are already stored.")

    for tip, resid in in_res.items():
        in_res[tip] = (resid - rmean) / rstdev

def edit_py_node(phy_in):
    """edit node to hep pep8 compliant"""
    for node in phy_in.get_terminals():
        node.comment = None

    for node in phy_in.get_nonterminals():
        if node.name is None and node.confidence:
            node.name = node.confidence
            node.confidence = None
        node.comment = None

def retrieve_genomes(
        by_lineage_in,
        known_seqs,
        ref_file,
        outgroup=None,
        earliest=True,
        callback=None):
    """
    Identify most recent sampled genome sequence for each Pangolin lineage.
    Export as FASTA for TreeTime analysis, including reference genome.

    :param by_lineage_in:  dict, return value from gisaid_utils::sort_by_lineage
    :param known_seqs:  dict, sequences used in Pango lineage designations
    :param ref_file:  str, path to FASTA file containing reference genome
    :param outgroup:  str, path to FASTA containing outgroup sequence (optional).
                      Defaults to None, where reference is used as outgroup.
    :param earliest:  bool, if False then use most recent genome as lineage representative
    :param callback:  optional, callback function

    :return:  dict, aligned genome keyed by header "|{lineage}|{coldate}"
    """
    # load and parse reference genome
    with open(ref_file, encoding='utf-8') as ref_handle:
        refseq = convert_fasta(ref_handle)[0][1]
        seqs = [refseq] # used by fasttree to root the tree

    # optionally load outgroup sequence
    if outgroup:
        with open(outgroup, encoding='utf-8') as ref_handle:
            seqs = [convert_fasta(ref_handle)[0][1]]

    # allocate lists
    coldates = [None]
    out_lineage = ['reference']

    # retrieve unaligned genomes from database
    for lin, records in by_lineage_in.items():
        records = covizu.utils.batch_utils.unpack_records(records)

        # filter records for lineage-defining genomes
        curated = filter(
            lambda r: r['covv_virus_name'].replace(
                'hCoV-19/',
                '').replace(
                ' ',
                '_')  # issue #313
            in known_seqs,
            records)

        curated = list(curated)  # resolve generator
        if len(curated) == 0:
            if callback:
                callback(
                    f"Error in retrieve_genomes(): no sequence names for lineage {lin} "
                    "in designated list; may need to update data/pango-designation/lineages.csv",
                    level='WARN')
            curated = records

        intermed = [(r['covv_collection_date'], r['diffs'], r['missing'])
                    for r in curated]

        parse_intermed(intermed, coldates, earliest, refseq, seqs)

        # update lists
        out_lineage.append(lin)

        # reconstruct aligned sequence from feature vector

    # generate new headers in {name}|{accession}|{date} format expected by
    # treetime()
    return dict(zip(map(lambda xy: f'|{xy[0]}|{xy[1]}', zip(out_lineage, coldates)), seqs))

def parse_intermed(intermed, coldates, earliest, refseq, seqs):
    """make retrieve genomes pep8 compliant"""
    intermed.sort(reverse=True)  # descending order
    coldate, diffs, missing = intermed[-1] if earliest else intermed[0]

    # update lists
    coldates.append(coldate)

    # reconstruct aligned sequence from feature vector
    seq = apply_features(diffs, missing=missing, refseq=refseq)
    seqs.append(seq)

def parse_args():
    """parse args"""
    parser = argparse.ArgumentParser(
        description="Generate inputs for TreeTime analysis."
    )
    parser.add_argument(
        'json',
        type=str,
        help='input, JSON produced by gisaid_utils.py')

    # issue #230
    path = covizu.__path__[0]
    if '.egg' in path:

        ref_file = pkg_resources.resource_filename(
            'covizu', 'data/NC_045512.fa')
    else:
        ref_file = os.path.join(path, "data", "NC_045512.fa")

    parser.add_argument('--ref', type=str, default=ref_file,
                        help="input, FASTA file with reference genome")
    parser.add_argument(
        '--misstol',
        type=int,
        default=300,
        help="optional, maximum tolerated number of missing bases per "
        "genome (default 300).")
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='optional, specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')
    parser.add_argument(
        '--earliest',
        action='store_true',
        help="option, select earliest genome per lineage; otherwise"
        " default to most recent samples.")

    parser.add_argument(
        '--outdir',
        default='data/',
        help='optional, directory to write TreeTime output files')
    parser.add_argument('--ft2bin', default='fasttree2',
                        help='optional, path to fasttree2 binary executable')
    parser.add_argument('--ttbin', default='treetime',
                        help='optional, path to treetime binary executable')
    parser.add_argument(
        '--lineages',
        type=str,
        default=os.path.join(
            covizu.__path__[0],
            "data/pango-designation/lineages.csv"),
        help="optional, path to CSV file containing Pango lineage designations.")

    parser.add_argument('--outfile', default='data/timetree.nwk',
                        help='output, path to write Newick tree string')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cb = Callback()

    cb.callback("Retrieving genomes")
    with open(args.json, encoding='utf-8') as handle:
        by_lineage = json.load(handle)

    cb.callback("Parsing Pango lineage designations")
    with open(args.lineages, encoding='utf-8') as handle:
        header = next(handle)
        if header != 'taxon,lineage\n':
            cb.callback(
                f"Error: {args.lineages} does not contain expected header row 'taxon,lineage'")
            sys.exit()
        lineages = {}
        for line in handle:
            try:
                taxon, lineage = line.strip().split(',')
                if taxon and lineage:
                    lineages.update({taxon: lineage})
                else:
                    cb.callback(
                        f"Warning '{line}': taxon or lineage is missing",
                        level='WARN')
            except ValueError as e:
                cb.callback(
                    f"Warning: There is an issue with the line '{line}' in lineages.csv, exit {e}",
                    level='WARN')

        cb.callback("Identifying lineage representative genomes")
        fasta = retrieve_genomes(
            by_lineage,
            known_seqs=lineages,
            ref_file=args.ref,
            earliest=args.earliest,
            callback=cb.callback)

        cb.callback(f"Reconstructing tree with {args.ft2bin}")
        nwk = fasttree(fasta, binpath=args.ft2bin)

        cb.callback(f"Reconstructing time-scaled tree with {args.ttbin}")
        nexus_file = treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                            clock=args.clock)

        timetree, residuals = parse_nexus(nexus_file, fasta)
        Phylo.write(timetree, file=args.outfile, format='newick')

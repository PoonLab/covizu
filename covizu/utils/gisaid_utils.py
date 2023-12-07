import lzma
import json
from datetime import date
import argparse
import os
import sys
import subprocess
from datetime import datetime
import getpass

import covizu
from covizu import minimap2
#from covizu.minimap2 import minimap2, encode_diffs
from covizu.utils.seq_utils import *
from covizu.utils.progress_utils import Callback

import gc


def download_feed(url, user, password):
    """
    Download xz file from GISAID.  Note this requires confidential URL, user and password
    information that we are not distributing with the source code.
    :param url:  str, address to retrieve xz-compressed provisioning file
    :param user:  str, GISAID username
    :param password:  str, access credentials - if None, query user
    :return:  str, path to time-stamped download file
    """
    if url is None:
        print("Error: no URL specified in download_feed()")
        sys.exit()
    if user is None:
        user = getpass.getpass("GISAID username: ")
    if password is None:
        password = getpass.getpass()
    timestamp = datetime.now().isoformat().split('.')[0]
    outfile = "data/provision.{}.json.xz".format(timestamp)
    subprocess.check_call(["wget", "--user", user, "--password", password, "-O", outfile, url])
    return outfile


def load_gisaid(path, minlen=29000, mindate='2019-12-01', callback=None,
                fields=("covv_accession_id", "covv_virus_name", "covv_lineage",
                        "covv_collection_date", "covv_location", "sequence"),
                debug=None
):
    """
    Read in GISAID feed as xz compressed JSON, applying some basic filters

    :param path:  str, path to xz-compressed JSON
    :param minlen:  int, minimum genome length
    :param mindate:  datetime.date, earliest reasonable sample collection date
    :param callback:  function, optional callback function
    :param fields:  tuple, fieldnames to keep
    :param debug:  int, if >0 then limits input JSON for debugging

    :yield:  dict, contents of each GISAID record
    """
    mindate = fromisoformat(mindate)
    rejects = {'short': 0, 'baddate': 0, 'nonhuman': 0, 'nolineage': 0}
    with lzma.open(path, 'rb') as handle:
        for ln, line in enumerate(handle):
            if debug and ln > debug:
                break
            record = json.loads(line)

            # remove unused data
            record = dict([(k, record[k]) for k in fields])

            qname = record['covv_virus_name'].strip().replace(',', '_').replace('|', '_')  # issue #206,#464
            country = qname.split('/')[1]
            if country == '' or country[0].islower():
                # reject mangled labels and non-human isolates
                rejects['nonhuman'] += 1
                continue

            record['covv_virus_name'] = qname  # in case we removed whitespace

            seq = record['sequence'].replace('\n', '')
            if len(seq) < minlen:
                # reject sequences that are too short
                rejects['short'] += 1
                continue
            record['sequence'] = seq

            if record['covv_collection_date'].count('-') != 2:
                # reject sequences without complete collection date
                rejects['baddate'] += 1
                continue
            coldate = fromisoformat(record['covv_collection_date'])
            if coldate is None or coldate < mindate or coldate > date.today():
                # reject sequences with nonsense collection date
                rejects['baddate'] += 1
                continue

            yield record

    if callback:
        callback("Rejected {short} short genomes\n"
                 "         {nolineage} with no lineage assignment\n"
                 "         {baddate} records with bad dates\n"
                 "         {nonhuman} non-human genomes".format(**rejects))


def batch_fasta(gen, cur, size=100):
    """
    Concatenate sequence records in stream into FASTA-formatted text in batches of
    <size> records.
    :param gen:  generator, return value of load_gisaid()
    :param size:  int, number of records per batch
    :yield:  str, list; FASTA-format string and list of records (dict) in batch
    """
    stdin = ''
    batch = []

    for i, record in enumerate(gen, 1):
        qname = record['covv_virus_name']
        sequence = record.pop('sequence')
        
        cur.execute("SELECT * FROM SEQUENCES WHERE qname = '%s'"%qname)
        result = cur.fetchone()
        if result:
            record.update({'diffs': result["diffs"], 'missing': result["missing"]})
        else:
            stdin += '>{}\n{}\n'.format(qname, sequence)
        batch.append(record)
        if i > 0 and i % size == 0:
            yield stdin, batch
            stdin = ''
            batch = []

    if batch:
        yield stdin, batch


def extract_features(batcher, ref_file, cur, binpath='minimap2', nthread=3, minlen=29000):
    """
    Stream output from JSON.xz file via load_gisaid() into minimap2
    via subprocess.

    :param batcher:  generator, returned by batch_fasta()
    :param ref_file:  str, path to reference genome (FASTA format)
    :param binpath:  str, path to minimap2 binary executable
    :param nthread:  int, number of threads to run minimap2
    :param minlen:  int, minimum genome length

    :yield:  dict, record augmented with genetic differences and missing sites;
    """
    with open(ref_file) as handle:
        reflen = len(convert_fasta(handle)[0][1])

    for fasta, batch in batcher:
        new_records = {}
        for record in batch:
            if 'diffs' in record:
                yield record
            else:
                new_records[record['covv_virus_name']] = record

        # If fasta is empty, no need to run minimap2
        if len(fasta) == 0:
            continue

        mm2 = minimap2.minimap2(fasta, ref_file, stream=True, path=binpath, nthread=nthread,
                       minlen=minlen)
        result = list(minimap2.encode_diffs(mm2, reflen=reflen))
        for qname, diffs, missing in result:
            # reconcile minimap2 output with GISAID record
            record = new_records[qname]
            record.update({'diffs': diffs, 'missing': missing})

            # inserting diffs and missing as json strings
            cur.execute("INSERT INTO SEQUENCES VALUES(%s, %s, %s, %s, %s, %s, %s)",
                        [json.dumps(v) if k in ['diffs', 'missing'] else v for k, v in record.items()])
            yield record


def filter_problematic(records, origin='2019-12-01', rate=0.0655, cutoff=0.005,
                       maxtime=1e3, vcf_file='data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf',
                       encoded=False,
                       misstol=300, callback=None):
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.

    :param records:  generator, records from extract_features() or encode_diffs()
    :param origin:  str, date of root sequence in ISO format (yyyy-mm-dd)
    :param rate:  float, molecular clock rate (subs/genome/day), defaults
                  to 8e-4 * 29900 / 365
    :param cutoff:  float, use 1-cutoff to compute quantile of Poisson
                    distribution, defaults to 0.005
    :param maxtime:  int, maximum number of days to cache Poisson quantiles
    :param vcf_file:  str, path to VCF file
    :param encoded: bool, flag to indicate if records are from encode_diffs() 
    :param misstol:  int, maximum tolerated number of uncalled bases
    :param callback:  function, option to print messages to console
    :yield:  generator, revised records
    """
    # load resources
    mask = load_vcf(vcf_file)
    qp = QPois(quantile=1-cutoff, rate=rate, maxtime=maxtime, origin=origin)

    n_sites = 0
    n_outlier = 0
    n_ambig = 0
    for record in records:
        if type(record) is not dict:
            qname, diffs, missing = record  # unpack tuple
        else:
            try:
                # loading json strings of old records
                record['diffs'] = json.loads(record['diffs'])
                record['missing'] = json.loads(record['missing'])
            except TypeError:
                # passing for new records
                pass
            finally:
                diffs = record['diffs']

        # exclude problematic sites
        filtered = []
        for typ, pos, alt in diffs:
            if typ == '~' and int(pos) in mask and alt in mask[pos]['alt']:
                continue
            if typ != '-' and 'N' in alt:
                # drop substitutions and insertions with uncalled bases
                continue
            filtered.append(tuple([typ, pos, alt]))

        ndiffs = len(filtered)
        n_sites += len(diffs) - ndiffs

        if not encoded:
            record['diffs'] = filtered

            # exclude genomes with excessive divergence from reference
            coldate = record['covv_collection_date']
            if qp.is_outlier(coldate, ndiffs):
                n_outlier += 1
                continue

            # exclude genomes with too much missing data
            if total_missing(record) > misstol:
                n_ambig += 1
                continue
            yield record
        else:
            if type(record) is dict:
                qname, missing = record['qname'], record['missing']

            yield [qname, filtered, missing]

    if callback:
        callback("filtered {} problematic features".format(n_sites))
        callback("         {} genomes with excess missing sites".format(n_ambig))
        callback("         {} genomes with excess divergence".format(n_outlier))


def sort_by_lineage(records, callback=None, interval=10000):
    """
    Resolve stream into a dictionary keyed by Pangolin lineage.
    Note: records yielded from generator accumulate in this function.

    :param records:  generator, return value of extract_features()
    :param callback:  optional, progress monitoring
    :param interval:  int, frequency to report alignment progress (genomes)
    :return:  dict, lists of records keyed by lineage
    """
    result = {}

    for i, record in enumerate(records):
        if callback and i % interval == 0:
            callback('aligned {} records'.format(i))

        lineage = record['covv_lineage']
        diffs = record.pop('diffs')  # REMOVE entry from record!
        if diffs is not None:
            diffs.sort()
        key = ','.join(['|'.join(map(str, diff)) for diff in diffs])

        if str(lineage) == "None" or lineage == '':
            # discard uncategorized genomes, #324, #335
            continue

        if lineage not in result:
            result.update({lineage: {}})
        if key not in result[lineage]:
            result[lineage].update({key: []})

        result[lineage][key].append(record)

    return result


def convert_json(infile, provision):
    """
    Convert clusters JSON file from main branch format to EpiCov format.
    :param infile:  open file stream ('r') to clusters JSON
    :param provision:  str, path to GISAID provision file
    :return:  dict, revised clusters to serialize to new JSON file
    """
    # first generate dictionary from provision keyed by accession
    metadata = {}
    with lzma.open(provision, 'rb') as handle:
        for line in handle:
            record = json.loads(line)
            metadata.update({record['covv_accession_id']: {
                'name': record['covv_virus_name'],
                'location': record['covv_location'],
                'coldate': record['covv_collection_date'],
                'gender': record['covv_gender'],
                'age': record['covv_patient_age'],
                'status': record['covv_patient_status']
            }})

    clusters = json.load(infile)
    for cluster in clusters:
        for variant, samples in cluster['nodes'].items():
            revised = []
            for coldate, accn, location, name in samples:
                md = metadata.get(accn, None)
                if md is None:
                    print("Failed to retrieve metadata for accession {}".format(accn))
                    sys.exit()
                revised.append([name, accn, location, coldate, md['gender'], md['age'], md['status']])

            # replace list of samples
            cluster['nodes'][variant] = revised

    return clusters


def parse_args():
    """ Command line help text"""
    parser = argparse.ArgumentParser("")
    parser.add_argument('outfile', type=argparse.FileType('w'), help="output, path to write JSON")

    parser.add_argument('--infile', type=str, default=None,
                        help="input, path to xz-compressed JSON")
    parser.add_argument('--url', type=str, 
                        help="URL to download provision file, defaults to environment variable.")
    parser.add_argument('--user', type=str, 
                        help="GISAID username, defaults to environment variable.")
    parser.add_argument('--password', type=str, 
                        help="GISAID password, defaults to environment variable.")

    parser.add_argument('--minlen', type=int, default=29000, help='option, minimum genome length')
    parser.add_argument('--mindate', type=str, default='2019-12-01',
                        help='option, earliest possible sample collection date (ISO format, default '
                             '2019-12-01)')

    parser.add_argument('--batchsize', type=int, default=500,
                        help='option, number of records to batch process with minimap2; limited by '
                             'buffer size for stdin redirection')

    parser.add_argument('--ref', type=str, help="option, path to reference genome (FASTA)",
                        default=os.path.join(covizu.__path__[0], "data/NC_045512.fa"))
    parser.add_argument('--binpath', type=str, default='minimap2',
                        help="option, path to minimap2 binary executable file")
    parser.add_argument('--mmthreads', type=int, default=8,
                        help='option, number of threads to run minimap2. Defaults to 8.')

    parser.add_argument("--vcf_file", type=str,
                        default=os.path.join(covizu.__path__[0], "data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf"),
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument("--debug", type=int, help="int, limit number of rows of input xz file to parse for debugging")

    args = parser.parse_args()

    if args.url is None and "GISAID_URL" in os.environ:
        args.url = os.environ["GISAID_URL"]
    if args.user is None and "GISAID_USER" in os.environ:
        args.user = os.environ["GISAID_USER"]
        # otherwise download_feed() will prompt for username
    if args.password is None and "GISAID_PSWD" in os.environ:
        args.password = os.environ["GISAID_PSWD"]    
        # otherwise download_feed() will prompt for password
    
    return args 


if __name__ == '__main__':
    args = parse_args()
    cb = Callback()

    cb.callback("Processing GISAID feed data")

    # download xz file if not specified by user
    if args.infile is None:
        args.infile = download_feed(args.url, args.user, args.password)

    loader = load_gisaid(args.infile, minlen=args.minlen, mindate=args.mindate,
                         debug=args.debug)
    batcher = batch_fasta(loader, size=args.batchsize)
    aligned = extract_features(batcher, ref_file=args.ref, binpath=args.binpath,
                               nthread=args.mmthreads, minlen=args.minlen)
    filtered = filter_problematic(aligned, vcf_file=args.vcf_file, callback=cb.callback)
    by_lineage = sort_by_lineage(filtered, callback=cb.callback)

    # serialize to JSON file
    json.dump(by_lineage, args.outfile)

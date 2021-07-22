from datetime import date
import bisect
import pkg_resources

from scipy.stats import poisson
from scipy.optimize import root


def read_seq(handle):
    """
    Read sequence from plain text file (no format).  Used for importing
    reference sequence.

    :param handle:
    :return:  str, sequence
    """
    seq = ''
    for line in handle:
        seq += line.strip()
    return seq


def parse_label(label):
    """
    Extract country and date of sample collection from GISAID label.
    Return date as None if collection date is ambiguous.

    :param label:  str, GISAID sequence label
    :return: (country, date)
    """
    info, epi_id, ymd = label.split('|')
    country = info.split('/')[1]
    try:
        year, month, day = list(map(int, ymd.split('-')))
        return country, date(year, month, day)
    except ValueError:
        return country, None
    except:
        raise


def iter_fasta(handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.

    :param handle:  open stream to FASTA file in read mode
    :yield tuples, (header, sequence)
    """
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    yield h, sequence


def convert_fasta(handle):
    """
    Parse FASTA file as a list of header, sequence list objects
    :param handle:  open file stream
    :return:  List of [header, sequence] records
    """
    result = []
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    result.append([h, sequence])  # handle last entry
    return result


def total_missing(row):
    """ Calculate the total number of missing sites from closed-open interval annotations """
    res = 0
    if type(row) is dict:
        missing = row['missing']
    else:
        _, _, missing = row

    for left, right in missing:
        res += right-left
    return res


def apply_features(diffs, missing, refseq):
    """
    Reconstitute genome sequence from feature vector (genetic differences) and
    missing data vector.

    :param row:  list, entry from features list returned by import_json()
    :param refseq:  str, reference genome
    :return:  str, aligned genome
    """
    result = list(refseq)  # strings are not mutable

    # apply missing intervals
    for left, right in missing:
        for i in range(left, right):
            result[i] = 'N'

    # apply substitutions and deletions (skip insertions)
    for dtype, pos, diff in diffs:
        if dtype == '~':
            result[pos] = diff
        elif dtype == '-':
            for i in range(pos, pos+diff):
                result[i] = '-'

    return ''.join(result)


def fromisoformat(dt):
    """ Convert ISO date to Python datetime.date object to support Python <3.7 """
    year, month, day = map(int, dt.split('-'))
    return date(year, month, day)


class QPois:
    """
    Cache the quantile transition points for Poisson distribution for a given
    rate <L> and varying time <t>, s.t. \exp(-Lt)\sum_{i=0}^{k} (Lt)^i/i! = Q.
    """
    def __init__(self, quantile, rate, maxtime, origin='2019-12-01'):
        """
        :param quantile:  float, cut distribution at q and 1-q
        :param rate:  float, molecular clock rate (genome / day)
        :param maxtime:  int, maximum number of days to cache
        :param origin:  str, x-intercept of trend in ISO format
        """
        if quantile <= 0 or quantile >= 1:
            print("ERROR: QPois quantile argument must be on closed interval (0,1).")
        self.upperq = quantile if quantile > 0.5 else (1-quantile)
        self.lowerq = 1-self.upperq  # TODO: store timepoints for lower quantiles
        self.rate = rate
        self.maxtime = maxtime
        self.origin = fromisoformat(origin)

        self.timepoints_upper, self.timepoints_lower = self.compute_timepoints()

    def objfunc(self, x, k, q):
        """ Use root-finding to find transition point for Poisson CDF """
        return q - poisson.cdf(k=k, mu=self.rate*x)

    def compute_timepoints(self, maxk=100):
        """ Store transition points until time exceeds maxtime """
        timepoints_upper = []
        timepoints_lower = []
        tu = 0
        tl = 0
        for k in range(maxk):
            res_upper = root(self.objfunc, x0=tu, args=(k, self.upperq, ))
            res_lower = root(self.objfunc, x0=tl, args=(k, self.lowerq, ))
            if not res_upper.success:
                print("Error in QPois: failed to locate root, q={} k={} rate={}".format(
                    self.upperq, k, self.rate))
                break
            if not res_lower.success:
                print("Error in QPois: failed to locate root, q={} k={} rate={}".format(
                    self.lowerq, k, self.rate))
                break
            tu = res_upper.x[0]
            tl = res_lower.x[0]
            if tu > self.maxtime:
                break
            timepoints_upper.append(tu)
            if tl < 0:
                break
            if tl <= self.maxtime:
                timepoints_lower.append(tl)
        return timepoints_upper, timepoints_lower

    def lookup(self, time, timepoints):
        """ Retrieve quantile count, given time """
        return bisect.bisect(timepoints, time)

    def is_outlier(self, coldate, ndiffs):
        if type(coldate) is str:
            coldate = fromisoformat(coldate)
        dt = (coldate - self.origin).days
        qmax = self.lookup(dt, self.timepoints_upper)
        qmin = self.lookup(dt, self.timepoints_lower)
        if ndiffs > qmax or ndiffs <= qmin:
            return True
        return False


def filter_outliers(iter, origin='2019-12-01', rate=0.0655, cutoff=0.005, maxtime=1e3):
    """
    Exclude genomes that contain an excessive number of genetic differences
    from the reference, assuming that the mean number of differences increases
    linearly over time and that the variation around this mean follows a
    Poisson distribution.

    :param iter:  generator, returned by encode_diffs()
    :param origin:  str, date of root sequence in ISO format (yyyy-mm-dd)
    :param rate:  float, molecular clock rate (subs/genome/day), defaults
                  to 8e-4 * 29900 / 365
    :param cutoff:  float, use 1-cutoff to compute quantile of Poisson
                    distribution, defaults to 0.005
    :param maxtime:  int, maximum number of days to cache Poisson quantiles
    :yield:  tuples from generator that pass filter
    """
    qp = QPois(quantile=1-cutoff, rate=rate, maxtime=maxtime, origin=origin)
    for qname, diffs, missing in iter:
        coldate = qname.split('|')[-1]
        if coldate.count('-') != 2:
            continue
        ndiffs = len(diffs)
        if qp.is_outlier(coldate, ndiffs):
            # reject genome with too many differences given date
            continue
        yield qname, diffs, missing


def load_vcf(vcf_file="data/problematic_sites_sarsCov2.vcf"):
    """
    Load VCF of problematic sites curated by Nick Goldman lab
    NOTE: The curators of this VCF used MN908947.3, which is identical to NC_045512.
    *** It is very important to check that your reference is compatible! ***
    TODO: align user's reference to NC_045512 to generate custom coordinate system

    :param vcf_file:  str, path to VCF file
    :return:  dict, tuples keyed by reference coordinate
    """
    vcf = open(vcf_file)
    mask = {}
    for line in vcf.readlines():
        if line.startswith('#'):
            continue
        try:
            _, pos, _, ref, alt, _, filt, info = line.strip().split('\t')
        except ValueError:
            raise
        if filt == 'mask':
            mask.update({int(pos)-1: {  # convert to 0-index
                'ref': ref, 'alt': alt, 'info': info}
            })
    return mask


def filter_problematic_sites(obj, mask, callback=None):
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.
    FIXME: this duplicates some functionality of filter_problematic(), #290

    :param obj:  list, entries are (1) dicts returned by import_json or (2) tuples
    :param mask:  dict, problematic site index from load_vcf()
    :param vcf_file:  str, path to VCF file
    :return:
    """
    # apply filters to feature vectors
    count = 0
    result = []
    for row in obj:
        if type(row) is dict:
            qname, diffs, missing = row['qname'], row['diffs'], row['missing']
        else:
            qname, diffs, missing = row  # unpack tuple

        filtered = []
        for typ, pos, alt in diffs:
            if typ == '~' and int(pos) in mask and alt in mask[pos]['alt']:
                continue
            if typ != '-' and 'N' in alt:
                # drop substitutions and insertions with uncalled bases
                continue
            filtered.append(tuple([typ, pos, alt]))

        count += len(diffs) - len(filtered)
        result.append([qname, filtered, missing])

    if callback:
        callback('filtered {} problematic features'.format(count))
    return result


class SC2Locator:
    """
    Annotate features
    """
    def __init__(self, ref_file=pkg_resources.resource_filename('covizu', 'data/NC_045512.fa')):
        self.gcode = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            '---': '-', 'XXX': '?'
        }
        self.orfs = {
            'orf1a': (265, 13468),
            'orf1b': (13467, 21555),
            'S': (21562, 25384),
            'orf3a': (25392, 26220),
            'E': (26244, 26472),
            'M': (26522, 27191),
            'orf6': (27201, 27387),
            'orf7a': (27393, 27759),
            'orf7b': (27755, 27887),
            'orf8': (27893, 28259),
            'N': (28273, 29533),
            'orf10': (29557, 29674)
        }

        # load reference genome
        with open(ref_file) as handle:
            fasta = convert_fasta(handle)
            self.refseq = fasta[0][1]

    def parse_mutation(self, feat):
        """
        Map feature from reference nucleotide coordinate system to amino
        acid substitutions, if relevant.
        :param feat:  tuple (str, int, str); type, position and alternate state
        :return:  str, AA or indel string; or None if synonymous nucleotide substitution
        """
        # unpack feature
        typ, pos, alt = feat
        if typ == '~':
            # is substitution within a reading frame?
            this_orf = None
            this_left, this_right = None, None
            for orf, coords in self.orfs.items():
                left, right = coords
                if left <= pos < right:
                    this_orf = orf
                    this_left, this_right = left, right
                    break

            # does the mutation change an amino acid?
            if this_orf:
                # retrieve codons
                codon_left = 3 * ((pos-this_left)//3)
                codon_pos = (pos-this_left) % 3

                rcodon = self.refseq[this_left:this_right][codon_left:(codon_left+3)]
                ramino = self.gcode[rcodon]

                qcodon = list(rcodon)
                qcodon[codon_pos] = alt
                qcodon = ''.join(qcodon)
                qamino = self.gcode[qcodon]

                if ramino != qamino:
                    return 'aa:{}:{}{:0.0f}{}'.format(this_orf, ramino, 1+codon_left/3, qamino)

        elif typ == '+':
            return 'ins:{}:{}'.format(pos+1, len(alt))

        elif typ == '-':
            return 'del:{}:{}'.format(pos+1, alt)

        # otherwise
        return None


def batch_fasta(gen, size=100):
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
        qname = record['label']
        sequence = record.pop('sequence')
        stdin += '>{}\n{}\n'.format(qname, sequence)
        batch.append(record)
        if i > 0 and i % size == 0:
            yield stdin, batch
            stdin = ''
            batch = []

    if batch:
        yield stdin, batch


def filter_problematic(records, origin='2019-12-01', rate=0.0655, cutoff=0.005,
                       maxtime=1e3, vcf_file='data/problematic_sites_sarsCov2.vcf',
                       misstol=300, callback=None):
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.

    :param records:  generator, records from extract_features()
    :param origin:  str, date of root sequence in ISO format (yyyy-mm-dd)
    :param rate:  float, molecular clock rate (subs/genome/day), defaults
                  to 8e-4 * 29900 / 365
    :param cutoff:  float, use 1-cutoff to compute quantile of Poisson
                    distribution, defaults to 0.005
    :param maxtime:  int, maximum number of days to cache Poisson quantiles
    :param vcf_file:  str, path to VCF file
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
        # exclude problematic sites
        filtered = []
        diffs = record['diffs']
        for typ, pos, alt in diffs:
            if typ == '~' and int(pos) in mask and alt in mask[pos]['alt']:
                continue
            if typ != '-' and 'N' in alt:
                # drop substitutions and insertions with uncalled bases
                continue
            filtered.append(tuple([typ, pos, alt]))

        ndiffs = len(filtered)
        n_sites += len(diffs) - ndiffs
        record['diffs'] = filtered

        # exclude genomes with excessive divergence from reference
        coldate = record['coldate']
        if qp.is_outlier(coldate, ndiffs):
            n_outlier += 1
            continue

        # exclude genomes with too much missing data
        if total_missing(record) > misstol:
            n_ambig += 1
            continue

        yield record

    if callback:
        callback("filtered {} problematic features".format(n_sites))
        callback("         {} genomes with excess missing sites".format(n_ambig))
        callback("         {} genomes with excess divergence".format(n_outlier))

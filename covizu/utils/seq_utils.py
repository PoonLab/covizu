from datetime import date
import bisect

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


def apply_features(row, refseq):
    """
    Reconstitute genome sequence from feature vector (genetic differences) and
    missing data vector.

    :param row:  list, entry from features list returned by import_json()
    :param refseq:  str, reference genome
    :return:  str, aligned genome
    """
    result = list(refseq)  # strings are not mutable
    _, diffs, missing = row

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
        self.q = quantile
        self.rate = rate
        self.maxtime = maxtime
        self.origin = fromisoformat(origin)

        self.timepoints = self.compute_timepoints()

    def objfunc(self, x, k):
        """ Use root-finding to find transition point for Poisson CDF """
        return self.q - poisson.cdf(k=k, mu=self.rate*x)

    def compute_timepoints(self, maxk=100):
        """ Store transition points until time exceeds maxtime """
        timepoints = []
        t = 0
        for k in range(maxk):
            res = root(self.objfunc, x0=t, args=(k, ))
            if not res.success:
                print("Error in QPois: failed to locate root, q={} k={} rate={}".format(
                    self.q, k, self.rate))
                print(res)
                break
            t = res.x[0]
            if t > self.maxtime:
                break
            timepoints.append(t)
        return timepoints

    def lookup(self, time):
        """ Retrieve quantile count, given time """
        return bisect.bisect(self.timepoints, time)

    def is_outlier(self, coldate, ndiffs):
        if type(coldate) is str:
            coldate = fromisoformat(coldate)
        dt = (coldate - self.origin).days
        qmax = self.lookup(dt)
        if ndiffs > qmax:
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
        _, pos, _, ref, alt, _, filt, info = line.strip().split()
        if filt == 'mask':
            mask.update({int(pos)-1: {  # convert to 0-index
                'ref': ref, 'alt': alt, 'info': info}
            })
    return mask


def filter_problematic(obj, mask, callback=None):
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.

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
            qname, diffs, missing = row

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

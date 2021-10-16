import urllib.request
import tarfile
import sys
import lzma
import gzip
from io import BytesIO
from datetime import datetime
import argparse
from csv import DictReader


def progress(msg):
    timestamp = datetime.now().isoformat()
    print("[{}] {}".format(timestamp, msg))


def retrieve_data(url, bufsize=16384):
    """
    Download data and save in memory
    @param url:  str, URL to download all data from VirusSeq API
    @param bufsize:  int, number of bytes to read per iteration, defaults to 2^14
    @return: BytesIO object
    """
    resp = urllib.request.urlopen(url)
    if resp.status != 200:
        print("Response status code {}".format(resp.status))
        sys.exit()

    tmpfile = BytesIO()
    while True:
        buf = resp.read(bufsize)
        if not buf:
            break  # end of content
        tmpfile.write(buf)
    resp.close()
    return tmpfile


def get_filenames(tmpfile):
    """
    Retrieve FASTA and TSV file names from tar
    @return (str, str)
    """
    tmpfile.seek(0)
    tf = tarfile.open(fileobj=tmpfile, mode="r|gz")

    files = tf.getnames()
    fasta_file = [f for f in files if f.endswith('.fasta')]
    if len(fasta_file) == 0:
        progress("ERROR: Failed to locate FASTA file in tar file")
        sys.exit()

    tsv_file = [f for f in files if f.endswith('.tsv')]
    if len(fasta_file) == 0:
        progress("ERROR: Failed to locate TSV file in tar file")
        sys.exit()
    tf.close()
    return fasta_file[0], tsv_file[0]


def export_files(fasta_file, tsv_file, tmpfile):
    """
    @param fasta_file:  str, name of FASTA file in tar
    @param tsv_file:  str, name of TSV file in tar
    @param tmpfile:  BytesIO, tarfile contents in memory
    @return str, str:  FASTA xz and TSV gz file names
    """
    timestamp = datetime.now().isoformat().split('.')[0]
    tmpfile.seek(0)
    tf = tarfile.open(fileobj=tmpfile, mode="r:gz")
    fasta = tf.extractfile(fasta_file)
    xzfile = "virusseq.{}.fasta.gzip".format(timestamp)
    with gzip.open(xzfile, 'w') as handle:
        handle.write(fasta.read())

    progress("exporting metadata TSV to gzip-compressed file")
    tsv = tf.extractfile(tsv_file)
    gzfile = "virusseq.{}.metadata.tsv.gz".format(timestamp)
    with gzip.open(gzfile, 'wb') as handle:
        handle.write(tsv.read())

    return xzfile, gzfile


def process_virusseq(fasta_file, tsv_file, tmpfile, fieldnames, callback=None):
    """
    Parse VirusSeq sequence and metadata feeds

    :param fasta_file:  str, file name of
    """
    # open metadata feed
    tmpfile.seek(0)
    tf = tarfile.open(fileobj=tmpfile, mode="r:gz")
    handle = tf.extractfile(tsv_file)

    # check fieldnames
    reader = DictReader(handle, delimiter='\t')

    for field in fieldnames:
        if field not in reader.fieldnames:
            if callback:
                callback("Missing fieldname {} in metadata CSV".format(field), level='ERROR')
                callback(reader.fieldnames, level='ERROR')
            sys.exit()

    # parse feed
    metadata = {}
    for row in reader:
        metadata.update({row['strain']: {
            'accession': row['genbank_accession'],
            'coldate': row['date'],
            'region': None,
            'country': row['country'],
            'division': None,
            'lineage': row['pango_lineage']
        }})

    print(len(metadata))

    fasta = tf.extractfile(fasta_file)
    for line in handle:
        print(line)
        break


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Download data from VirusSeq portal database,"
                                     "compressing streams to reduce file size.")
    parser.add_argument("--url", default="https://singularity.virusseq-dataportal.ca/download/archive/all",
                        type=str, help="URL to retrieve all data via VirusSeq API")
    args = parser.parse_args()

    progress("Caching download in memory")
    tmpfile = retrieve_data(url=args.url)

    progress("Extracting file names in tar")
    fasta_file, tsv_file = get_filenames(tmpfile)

    progress("Exporting compressed files")
    paths = export_files(fasta_file, tsv_file, tmpfile)
    progress("Wrote outputs to {0} and {1}".format(*paths))

    sys.exit()
    # for VirusSeq feed
    fieldnames = ['isolate', 'specimen collector sample ID', 'sample collection date',
                  'geo_loc_name (country)', 'pango_lineage']
    process_virusseq(fasta_file, tsv_file, tmpfile, fieldnames)

    # for Nextstrain feed
    fieldnames = ['strain', 'genbank_accession', 'date', 'country', 'pango_lineage']
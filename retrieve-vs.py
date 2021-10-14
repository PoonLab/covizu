import urllib.request
import tarfile
import sys
import lzma
import gzip
from io import BytesIO

url = "https://singularity.virusseq-dataportal.ca/download/archive/all"

resp = urllib.request.urlopen(url)
if resp.status != 200:
    print("Response status code {}".format(resp.status))
    sys.exit()

# store download in temporary file
tmpfile = BytesIO()
while True:
    s = resp.read(16384)
    if not s:
        break
    tmpfile.write(s)
resp.close()

tmpfile.seek(0)
tf = tarfile.open(fileobj=tmpfile, mode="r|gz")

files = tf.getnames()
fasta_file = [f for f in files if f.endswith('.fasta')]
if len(fasta_file) == 0:
    print("ERROR: Failed to locate FASTA file in tar file")
    sys.exit()

tsv_file = [f for f in files if f.endswith('.tsv')]
if len(fasta_file) == 0:
    print("ERROR: Failed to locate TSV file in tar file")
    sys.exit()
tf.close()

tmpfile.seek(0)
tf = tarfile.open(fileobj=tmpfile, mode="r:gz")
fasta = tf.extractfile(fasta_file[0])
# TODO: timestamp the output filenames
with lzma.open("virusseq.fasta.xz", 'w') as handle:
    handle.write(fasta.read())

tsv = tf.extractfile(tsv_file[0])
with gzip.open("virusseq.metadata.tsv", 'wb') as handle:
    handle.write(tsv.read())

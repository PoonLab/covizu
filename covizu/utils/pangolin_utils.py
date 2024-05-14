import pandas as pd
import sklearn
import joblib

import argparse
import sys
import pkg_resources
import os
from csv import DictWriter

import covizu
from covizu.utils import seq_utils
from covizu.minimap2 import minimap2, stream_fasta


class Pangolin:
    def __init__(self, header_file, model_file):
        model_headers = joblib.load(header_file)
        self.indices = model_headers[1::]  # list of integers
        self.model = joblib.load(model_file)
        self.categories = ['A', 'C', 'G', 'T', 'N', '-']

    def classify(self, seq):
        """ Assign genome to one or more lineages """
        # convert sequence into list
        seqlist = [nt if nt in 'ACGT-' else 'N' for i,
                   nt in enumerate(seq) if i in self.indices]
        df = pd.DataFrame([seqlist], columns=self.indices)

        # add extra rows to ensure all categories are represented
        for nt in self.categories:
            df.loc[len(df)] = [nt] * len(self.indices)

        df = pd.get_dummies(df, columns=self.indices)  # one-hot encoding

        df.drop(df.tail(len(self.categories)).index,
                inplace=True)  # remove fake data

        predictions = self.model.predict_proba(df)
        return predictions

    def process_fasta(self, handle, ref_file, binpath, nthread, minlen=29000):
        """
        Run all genomes in local FASTA file through Pangolin
        :param handle:  file object, opened in read mode
        :param ref_file:  str, path to FASTA file containing reference genome
        :param binpath:  str, path to minimap2 binary
        :param nthread:  int, number of threads to run minimap2
        :param minlen:  int, reject genomes below this threshold
        :yield:  tuple, header and lineage
        """
        reflen = len(seq_utils.convert_fasta(open(ref_file))[0][1])
        mm2 = minimap2(
            handle,
            ref_file,
            stream=False,
            path=binpath,
            nthread=nthread,
            minlen=minlen)
        for header, aligned in stream_fasta(mm2, reflen=reflen):
            lineage = self.classify(aligned)
            yield header, lineage


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run Pangolin SARS-CoV-2 lineage classifier on FASTA inputs')

    parser.add_argument("fasta", type=argparse.FileType(
        "r"), help="input, path to FASTA file with genomes to process")
    parser.add_argument("outfile", type=argparse.FileType('w'),
                        help="output, path to write CSV output")

    # minimap2 arguments
    parser.add_argument('--mm2bin', type=str, default='minimap2',
                        help='str, path to minimap2 binary executable')
    parser.add_argument(
        '--ref',
        type=str,
        help="<input> path to target FASTA (reference)",
        default=os.path.join(
            covizu.__path__[0],
            "data/NC_045512.fa"))
    parser.add_argument('--minlen', help="<option> minimum sequence length, "
                                         "defaults to 29000nt.",
                        type=int, default=29000)
    parser.add_argument('-t', '--thread', type=int, default=3,
                        help="<option> number of threads")

    # allow user to specify custom file locations
    parser.add_argument(
        "--model_file",
        type=str,
        default=pkg_resources.resource_filename(
            'pangoLEARN',
            'data/decisionTree_v1.joblib'),
        help='str, path to PangoLEARN decisionTree_V1.joblib file')
    parser.add_argument(
        "--header_file",
        type=str,
        default=pkg_resources.resource_filename(
            'pangoLEARN',
            'data/decisionTreeHeaders_v1.joblib'),
        help='str, path to PangoLEARN decisionTreeHeaders_v1.joblib file')

    return parser.parse_args()


if __name__ == '__main__':
    """ Command-line interface """
    args = parse_args()
    pangolin = Pangolin(args.header_file, args.model_file)
    output = pangolin.process_fasta(
        args.fasta, ref_file=args.ref, binpath=args.mm2bin,
        nthread=args.thread, minlen=args.minlen
    )
    writer = DictWriter(args.outfile, fieldnames=['header', 'lineage'])
    for header, lineage in output:
        writer.writerow({'header': header, 'lineage': lineage})

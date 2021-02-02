"""
Modified script from PangoLEARN.py from PANGOLIN project
https://github.com/cov-lineages/pangolin
"""

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
import datetime

import joblib
import argparse
import pangoLEARN
import re
import sys

from covizu.utils.seq_utils import convert_fasta
from covizu.utils.db_utils import *


def parse_args():
    parser = argparse.ArgumentParser(description='pangorider. A derivative of PANGOLIN')
    site_packages = next(p for p in sys.path if 'site-packages' in p)
    parser.add_argument('--pangolindir', help ='PangoLEARN data dir, defaults to .../site_packages/panoLEARN/data/',
        default=site_packages+'/pangoLEARN/data/')
    parser.add_argument("--header_file", action="store", type=str, dest="header_file")
    parser.add_argument("--model_file", action="store", type=str, dest="model_file")
    parser.add_argument("--fasta", action="store", type=str, dest="sequences_file")
    parser.add_argument("-o","--outfile", action="store", type=str, dest="outfile")
    parser.add_argument('--indicies', help='Indicies to keep, defaults to [265:29674]')
    parser.add_argument('--filterout', default = 'debug/filtered.log', type =argparse.FileType('w+'),
        help='Log for filtered seqs')

    return parser.parse_args()


def build_sequence_handle(sequence_file):
    """
    Takes a cleaned and aligned fasta (sequence_file)
    """
    sequence_handle = []
    with open(sequence_file) as f:
        currentSeq = ""
        seqID = f.readline().strip()[1:]
        for line in f:
            if line.startswith('>'):
                #this is a header
                seqID = line.strip()[1:]
                sequence_handle.append([seqID, currentSeq])
            else:
                currentSeq+= line.strip()

    return sequence_handle


# function for handling weird sequence characters
def clean(x):
    if x == 'A' or x == 'C' or x == 'T' or x == '-':
        return x
    return 'N'


# generates data line
def encodeSeq(seq, indiciesToKeep):
    dataLine = []
    for i in indiciesToKeep:
        if i < len(seq):
            dataLine.extend(clean(seq[i]))

    return dataLine


# reads in the two data files
def readInAndFormatData(sequence_handle, indiciesToKeep, blockSize=1000):
    """
    sequence_handle will be a python list object [[header, seq], ... , [header, seq]]
    """
    idList = []
    seqList = []

    currentSeq = ""

    # parse the handle collecting a list of the sequences
    for seqid, currentSeq in sequence_handle:        
        idList.append(seqid)
        seqList.append(encodeSeq(currentSeq, indiciesToKeep))
        currentSeq = ""

        if len(seqList) == blockSize:
            yield idList, seqList
            idList = []
            seqList = []

        # gotta get the last one
        idList.append(seqid)
        seqList.append(encodeSeq(currentSeq, indiciesToKeep))

    yield idList, seqList


def filter_seqs(sequence_handle, outfile, max_prop_n=0.05, minlen=29000):
    """
    Filter FASTA file for partial and non-human SARS-COV-2 genome sequences.

    *Adapted from filtering.py*

    :param fasta_file:  open file stream to GISAID FASTA file
    :param outfile:  open file stream to write filtered FASTA file
    :param trim_left:  int, number of bases to drop from left
    :param trim_right:  int, number of bases to drop from right
    :param max_prop_n:  float, maximum proportion of N's (ambiguous bases)
                        tolerated per genome
    :param minlen:  int, minimum tolerated sequence length

    :return:  dict, containing lists of headers for rejected genomes
    """
    # lower-case label in place of country identifies non-human samples
    pat = re.compile('^[^/]+/[a-z]')
    pat2 = re.compile("^[HhCcOoVv]+-19/[A-Z][^/]+/[^/]+/[0-9-]+\|[^|]+\|[0-9]{4}-[0-9]+-[0-9]+")
    pat3 = re.compile('^-*')
    pat4 = re.compile('-*$')

    accessions = {}
    discards = {'nonhuman': [], 'ambiguous': [], 'short': [],
                'duplicates': [], 'mangled header': []}
    filter_handle = []

    # Set date vars, headers must be within this, otherwise consider it mangled
    today = datetime.datetime.today().date()
    init_day = datetime.date(2019, 12, 1)  # earliest possible Covid seq cannot be before Dec 2019

    for h, s in sequence_handle:
        if not type(h) == str or not type(s) == str:
            print("Error: entry {} not string type: sequence {}".format(h, s))
            continue
        if pat.findall(h):
            discards['nonhuman'].append(h)
            continue

        if len(s) < minlen:
            discards['short'].append(h)
            continue

        # apply sequence trims :Depreciated:
        #seq = s[trim_left:(-trim_right)]
        seq = s

        # this is conservative - all internal gaps are interpreted as deletions
        gap_prefix = len(pat3.findall(seq)[0])
        gap_suffix = len(pat4.findall(seq)[0])
        seqlen = len(seq) - gap_prefix - gap_suffix

        n_ambig = seq.count('?') + seq.count('N') + gap_prefix + gap_suffix
        if n_ambig / float(len(seq)) > max_prop_n:
            discards['ambiguous'].append(h)
            continue

        if pat2.search(h) is None:
            discards['mangled header'].append(h)
            continue

        desc, accn, coldate = h.split('|')
        if accn in accessions:
            discards['duplicates'].append(h)
            continue
        accessions.update({accn: desc})

        #check if date in proper format
        coldate_parsed = coldate.split('-')
        if len(coldate_parsed) != 3:
            discards['mangled header'].append(h)
            continue
        else:
            coldate_parsed = [int(i) for i in coldate_parsed]
            seq_date = datetime.date(coldate_parsed[0], coldate_parsed[1], coldate_parsed[2])
            if seq_date < init_day or seq_date > today:
                discards['mangled header'].append(h)
                continue

        # add genome to filtered handle
        filter_handle.append([h,s])

    outfile.write("Discarded {} non-human sequences.\n".format(len(discards['nonhuman'])))
    for h in discards['nonhuman']:
        outfile.write('  {}\n'.format(h))

    outfile.write("Discarded {} problematic sequences with prop N > {}.\n".format(
        len(discards['ambiguous']), max_prop_n
    ))
    for h in discards['ambiguous']:
        outfile.write('  {}\n'.format(h))

    outfile.write("Discarded {} sequences of length < {}\n".format(
        len(discards['short']), minlen
    ))
    for h in discards['short']:
        outfile.write('  {}\n'.format(h))

    outfile.write("Discarded {} sequences with duplicate accession numbers\n".format(
        len(discards['duplicates'])
    ))
    for h in discards['duplicates']:
        outfile.write('  {}\n'.format(h))

    outfile.write("Discarded {} sequences with mangled headers / ambiguous sample dates\n".format(
        len(discards['mangled header'])
    ))
    for h in discards['mangled header']:
        outfile.write('  {}\n'.format(h))

    return filter_handle


def classify_and_insert(header_file, model_file, sequence_handle, indiciesToKeep, database):
    """
    Re-written wrapper for pangoLEARN
    Assumes raw seqs in alignedseq are filtered (passed_qc) and aligned to LR757995.1
    Assumes that sequences in sequence_handle python list are tuples [header, alignedseq]

    """
    # loading the list of headers the model needs.
    model_headers = joblib.load(header_file)
    indiciesToKeep = model_headers[1:]

    print("loading model " + datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
    loaded_model = joblib.load(model_file)

    # write predictions to the LINEAGE table
    cursor, conn = open_connection(database)
    for idList, seqList in readInAndFormatData(sequence_handle, indiciesToKeep):
        print("processing block of {} sequences {}".format(
            len(seqList), datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
        ))

        # create a data from from the seqList
        df = pd.DataFrame(seqList, columns=indiciesToKeep)

        # possible nucleotide symbols
        categories = ['A', 'C', 'G', 'T', 'N', '-']

        # add extra rows to ensure all of the categories are represented, as otherwise 
        # not enough columns will be created when we call get_dummies
        for i in categories:
            line = [i] * len(indiciesToKeep)
            df.loc[len(df)] = line

        # get one-hot encoding
        df = pd.get_dummies(df, columns=indiciesToKeep)

        # get rid of the fake data we just added
        df.drop(df.tail(len(categories)).index,inplace=True)

        predictions = loaded_model.predict_proba(df)

        for index in range(len(predictions)):
            maxScore = 0
            maxIndex = -1

            # get the max probability score and its assosciated index
            for i in range(len(predictions[index])):
                if predictions[index][i] > maxScore:
                    maxScore = predictions[index][i]
                    maxIndex = i

            score = maxScore
            prediction = loaded_model.classes_[maxIndex]
            accession = idList[index].split('|')[1]
            payload = [accession, prediction, score, pangoLEARN.__version__, 'passed_qc', '']
            insert_lineage(cursor, payload)

    conn.commit()


if __name__ == '__main__':
    args = parse_args()
    sequence_handle = convert_fasta(open(args.sequences_file, 'r'))

    filtered_handle = filter_seqs(sequence_handle, args.filterout, max_prop_n=0.05, minlen=29000)

    if args.header_file is None and args.model_file is None:
        classify_and_insert(args.pangolindir + 'decisionTreeHeaders_v1.joblib',
                            args.pangolindir + 'decisionTree_v1.joblib', filtered_handle,
                            args.indicies, args.db)
    else:
        classify_and_insert(args.header_file, args.model_file, filtered_handle, args.indicies, args.db)


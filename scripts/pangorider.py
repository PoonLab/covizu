#Modified script from PangoLEARN.py from PANGOLIN project https://github.com/cov-lineages/pangolin
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
from datetime import datetime
import argparse
import joblib
import argparse
import pangoLEARN




def parse_args():
    parser = argparse.ArgumentParser(description='pangoLEARN.')
    parser.add_argument("--header-file", action="store", type=str, dest="header_file")
    parser.add_argument("--model-file", action="store", type=str, dest="model_file")
    parser.add_argument("--fasta", action="store", type=str, dest="sequences_file")
    parser.add_argument("-o","--outfile", action="store", type=str, dest="outfile")
    return parser.parse_args()


#Functions for converting raw sequences to acceptable handle for pangorider_main()
def filter_raw():
    """
    :TODO:
    Quality control step from pangolin combined with covizu filters 
    """
    pass

def align_raw():
    """
    :TODO:
    Need to implement minimap 
    Returns python list with 
    """
    pass
    #return handle

def build_sequence_handle(sequence_file):
    """
    Takes a cleaned and aligned fasta (sequence_file)
    """
    sequence_handle = []
    with open(sequencesFile) as f:
        currentSeq = ""
        seqID= f.readline().strip()[1:]
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


def classify_and_insert(header_file, model_file, sequence_handle, database):
    """
    Re-written wrapper for pangoLEARN
    Assumes raw seqs in alignedseq are filtered (passed_qc) and aligned to LR757995.1
    Assumes that sequences in sequence_handle python list are tuples [header, alignedseq]

    """
    # loading the list of headers the model needs.
    model_headers = joblib.load(header_file)
    indiciesToKeep = model_headers[1:]

    print("loading model " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
    loaded_model = joblib.load(model_file)

    # write predictions to the LINEAGE table
    cursor, conn = open_connection(database)
    for idList, seqList in readInAndFormatData(sequence_handle, indiciesToKeep):
        print("processing block of {} sequences {}".format(
            len(seqList), datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
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

if __name_ == '__main__':

    args = parse_args()
    sequence_handle = open(args.sequences_file, 'r')
    pangolin_version =
    pangorider_main(args.header_file, args.model_file, args.sequence_handle)
# copied directly from ArtPoon/gotoh2/gotoh_utils.py
import re

# regex for terminal gaps
tgap_regex = re.compile("^(-*)[^-].+[^-](-*)$")

# conversion from codon triplets to amino acid symbols
codon_dict = {
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

mixture_regex = re.compile('[WRKYSMBDHVN-]')

mixture_dict = {
    'W': 'AT', 'R': 'AG', 'K': 'GT', 'Y': 'CT', 'S': 'CG',
    'M': 'AC', 'V': 'AGC', 'H': 'ATC', 'D': 'ATG',
    'B': 'TGC', 'N': 'ATGC', '-': 'ATGC'
}

ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.items())


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
        if line.startswith('$'):  # skip comment
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()

    result.append([h, sequence])  # handle last entry
    return result


def len_terminal_gaps(seq):
    """
    Calculate lengths of terminal (prefix, suffix) gaps
    :param seq:  str, gapped nucleotide/amino acid sequence
    :return:  tuple, lengths of gap prefix and suffix
    """
    m = tgap_regex.findall(seq)
    return tuple(map(len, m[0]))


def get_boundaries(seq):
    """
    Legacy code.
    :param seq:  str, sequence to check for terminal gaps
    :return:  tuple, indices to slice seq to remove terminal gaps
    """
    left, right = len_terminal_gaps(seq)
    return left, len(seq) - right


def translate_nuc(seq, offset, resolve=False, return_list=False):
    """
    Translate nucleotide sequence into amino acid sequence.
    offset by X shifts sequence to the right by X bases
    Synonymous nucleotide mixtures are resolved to the corresponding residue.
    Nonsynonymous nucleotide mixtures are encoded with '?'

    :param seq:  string, nucleotide sequence
    :param offset:  integer, number of gaps to pad sequence on the left (5' end)
                    before translating
    :param resolve:  if True, attempt to convert ambiguous codons into an amino acid
    :param return_list:  if True, ambiguous codons that resolve to two or more amino
                         acids are returned as list objects.
    :return:  str if return_list is False (default), else list
    """

    seq = '-' * offset + seq

    aa_list = []
    aa_seq = ''  # use to align against reference, for resolving indels

    # loop over codon sites in nucleotide sequence
    for codon_site in range(0, len(seq), 3):
        codon = seq[codon_site:codon_site + 3]

        if len(codon) < 3:
            break

        # note that we're willing to handle a single missing nucleotide as an ambiguity
        if codon.count('-') > 1 or '?' in codon:
            if codon == '---':  # don't bother to translate incomplete codons
                aa_seq += '-'
                aa_list.append(['-'])
            else:
                aa_seq += '?'
                aa_list.append(['?'])
            continue

        # look for nucleotide mixtures in codon, resolve to alternative codons if found
        num_mixtures = len(mixture_regex.findall(codon))

        if num_mixtures == 0:
            aa = codon_dict[codon]
            aa_seq += aa
            aa_list.append([aa])

        elif num_mixtures == 1:
            resolved_AAs = []
            for pos in range(3):
                if codon[pos] in mixture_dict.keys():
                    for r in mixture_dict[codon[pos]]:
                        rcodon = codon[0:pos] + r + codon[(pos + 1):]
                        if codon_dict[rcodon] not in resolved_AAs:
                            resolved_AAs.append(codon_dict[rcodon])

            aa_list.append(resolved_AAs)

            if len(resolved_AAs) > 1:
                if resolve:
                    # for purposes of aligning AA sequences
                    # it is better to have one of the resolutions
                    # than a completely ambiguous '?'
                    aa_seq += resolved_AAs[0]
                else:
                    aa_seq += '?'
            else:
                aa_seq += resolved_AAs[0]

        else:
            aa_seq += '?'
            aa_list.append(['?'])

    if return_list:
        return aa_list

    return aa_seq


def apply_prot_to_nuc(aligned_prot, nuc):
    """
    Use aligned protein sequence to update the corresponding nucleotide
    sequence.
    :param aligned_prot:  str, aligned amino acid sequence
    :param nuc:  str, original nucleotide sequence
    :return: str, nucleotide sequence with codon gaps
    """
    res = ''
    i = 0
    for aa in aligned_prot:
        if aa == '-':
            res += '---'
            continue
        res += nuc[i:(i + 3)]
        i += 3
    return res
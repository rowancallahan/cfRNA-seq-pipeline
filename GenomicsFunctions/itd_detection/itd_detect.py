#!/usr/bin/env python

# VERSION: 0.1.0

import sys
from Bio import SeqIO
from itertools import groupby
from operator import itemgetter

def parse_fasta(handle_fasta):
    """
    Return a SeqIO structure containing FASTA sequence.
    """
    ref_dict = {}
    for record in SeqIO.parse(handle_fasta, "fasta"):
        ref_dict[record.id] = record.seq

    return ref_dict


def parse_fastq(handle_fastq):
    """
    Get sequences from FASTQ input.
    """
#    fastq_dict = {}
    for record in SeqIO.parse(handle_fastq, "fastq-sanger"):
#        fastq_dict[record.id] = record.seq
#        if "TCTCTTGGAAACTCCCATTTGAGATCATA" in record.seq:
#        if "GCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTTGCAAACTCTAATTTTCTCTGGAACTCCCATTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTATCATTATCTGAGGAGCCGGTCACCTGT" in record.seq:
        yield record.seq

#    return fastq_dict


def chunk_seq(seq, chunk=17, step=3):
    """
    Return chunks of the FASTA template sequence.
    """
    for b in range(0, len(seq)-chunk, step):
        yield seq[b:b+17]


def find_matches(handle_fastq, template):
    """
    """
    for record in parse_fastq(handle_fastq):
        for chunk in chunk_seq(template):
            if record.seq.count(chunk) > 1:
                yield record.seq
                break

def max_size_dups(template):
    """
    """
    for i in range(len(template)/2, 15, -1):
        for j in range(len(template)-(len(template)/2)):
            if template.count(template[j:i+j]) > 1:
                return template[j:i+j]
            

def longest_common_substring(text):
    """Get the longest common substrings and their positions.
    >>> longest_common_substring('banana')
    {'ana': [1, 3]}
    >>> text = "not so Agamemnon, who spoke fiercely to "
    >>> sorted(longest_common_substring(text).items())
    [(' s', [3, 21]), ('no', [0, 13]), ('o ', [5, 20, 38])]

    This function can be easy modified for any criteria, e.g. for searching ten
    longest non overlapping repeated substrings.
    """
    sa, rsa, lcp = suffix_array(text)
    maxlen = max(lcp)
    result = {}
    for i in range(1, len(text)):
        if lcp[i] == maxlen:
            j1, j2, h = sa[i - 1], sa[i], lcp[i]
            assert text[j1:j1 + h] == text[j2:j2 + h]
            substring = text[j1:j1 + h]
            if not substring in result:
                result[substring] = [j1]
            result[substring].append(j2)
    return dict((k, sorted(v)) for k, v in result.items())

def suffix_array(text, _step=16):
    """Analyze all common strings in the text.

    Short substrings of the length _step a are first pre-sorted. The are the 
    results repeatedly merged so that the garanteed number of compared
    characters bytes is doubled in every iteration until all substrings are
    sorted exactly.

    Arguments:
        text:  The text to be analyzed.
        _step: Is only for optimization and testing. It is the optimal length
               of substrings used for initial pre-sorting. The bigger value is
               faster if there is enough memory. Memory requirements are
               approximately (estimate for 32 bit Python 3.3):
                   len(text) * (29 + (_size + 20 if _size > 2 else 0)) + 1MB

    Return value:      (tuple)
      (sa, rsa, lcp)
        sa:  Suffix array                  for i in range(1, size):
               assert text[sa[i-1]:] < text[sa[i]:]
        rsa: Reverse suffix array          for i in range(size):
               assert rsa[sa[i]] == i
        lcp: Longest common prefix         for i in range(1, size):
               assert text[sa[i-1]:sa[i-1]+lcp[i]] == text[sa[i]:sa[i]+lcp[i]]
               if sa[i-1] + lcp[i] < len(text):
                   assert text[sa[i-1] + lcp[i]] < text[sa[i] + lcp[i]]
    >>> suffix_array(text='banana')
    ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [0, 1, 3, 0, 0, 2])

    Explanation: 'a' < 'ana' < 'anana' < 'banana' < 'na' < 'nana'
    The Longest Common String is 'ana': lcp[2] == 3 == len('ana')
    It is between  tx[sa[1]:] == 'ana' < 'anana' == tx[sa[2]:]
    """
    tx = text
    size = len(tx)
    step = min(max(_step, 1), len(tx))
    sa = list(range(len(tx)))
    sa.sort(key=lambda i: tx[i:i + step])
    grpstart = size * [False] + [True]  # a boolean map for iteration speedup.
    # It helps to skip yet resolved values. The last value True is a sentinel.
    rsa = size * [None]
    stgrp, igrp = '', 0
    for i, pos in enumerate(sa):
        st = tx[pos:pos + step]
        if st != stgrp:
            grpstart[igrp] = (igrp < i - 1)
            stgrp = st
            igrp = i
        rsa[pos] = igrp
        sa[i] = pos
    grpstart[igrp] = (igrp < size - 1 or size == 0)
    while grpstart.index(True) < size:
        # assert step <= size
        nextgr = grpstart.index(True)
        while nextgr < size:
            igrp = nextgr
            nextgr = grpstart.index(True, igrp + 1)
            glist = []
            for ig in range(igrp, nextgr):
                pos = sa[ig]
                if rsa[pos] != igrp:
                    break
                newgr = rsa[pos + step] if pos + step < size else -1
                glist.append((newgr, pos))
            glist.sort()
            for ig, g in groupby(glist, key=itemgetter(0)):
                g = [x[1] for x in g]
                sa[igrp:igrp + len(g)] = g
                grpstart[igrp] = (len(g) > 1)
                for pos in g:
                    rsa[pos] = igrp
                igrp += len(g)
        step *= 2
    del grpstart
    # create LCP array
    lcp = size * [None]
    h = 0
    for i in range(size):
        if rsa[i] > 0:
            j = sa[rsa[i] - 1]
            while i != size - h and j != size - h and tx[i + h] == tx[j + h]:
                h += 1
            lcp[rsa[i]] = h
            if h > 0:
                h -= 1
    if size > 0:
        lcp[0] = 0
    return sa, rsa, lcp


def check_repeats(string, frac=0.75):

    for i in range(2,5):
        curr_seq = string[:i]
        thresh = int((len(string)/i)*frac)
        if string.count(curr_seq) >= thresh:
            return 0


def main():

    # http://stackoverflow.com/questions/13560037/effcient-way-to-find-longest-duplicate-string-for-python-from-programming-pearl
    # Original implementation horribly slow.  If I can find my original sa code, will replace with this.

    # Handle reverse complements also?
    test_temp = "GCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTTGCAAACTCTAATTTTCTCTGGAACTCCCATTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTATCATTATCTGAGGAGCCGGTCACCTGT"

    # Set in argparse.
    chunk = 17
    step = 3
    
    handle_fasta = open(sys.argv[1], 'rU')
    handle_fastq = open(sys.argv[2], 'rU')

    ref_dict = parse_fasta(handle_fasta)
#    fastq_dict = parse_fastq(handle_fastq)
#    flt_e14 = ref_dict["13:28608148-28608422"]
    # For whole FLT3 gene.
    for ref_seq in ref_dict.itervalues():
        for fq_seq in parse_fastq(handle_fastq):
            longest = longest_common_substring(str(fq_seq))
            for long_repeat in longest:
                if len(long_repeat) > chunk and len(longest) == 1:
                    if long_repeat in ref_seq and check_repeats(long_repeat) != 0:
                        print(long_repeat)
                else:
                    break
#    for output in find_matches(handle_fastq, flt_e14):
#        print(max_size_dups(output))


if __name__ == "__main__":
    main()


#print(longest_common_substring("GCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTTGCAAACTCTAATTTTCTCTGGAACTCCCATTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTATCATTATCTGAGGAGCCGGTCACCTGT"))




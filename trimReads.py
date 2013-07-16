#!/usr/bin/python

# Copyright (c) 2013, Stephen Fisher and Junhyong Kim, University of
# Pennsylvania.  All Rights Reserved.
#
# You may not use this file except in compliance with the Kim Lab License
# located at
#
#     http://kim.bio.upenn.edu/software/LICENSE
#
# Unless required by applicable law or agreed to in writing, this
# software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied.  See the License
# for the specific language governing permissions and limitations
# under the License.

"""
by: S. Fisher, 2013

usage: trimReads.py [-h] [-v] [-p] [-m MIN_LEN] [-c3 NUM_CUT_3] [-c5 NUM_CUT_5] [-rN] [-rAT REMOVE_AT] [-c CONTAMINANTS_FA] -f
                    FORWARD_FQ [-r REVERSE_FQ] -o OUTPUT_PREFIX

Contaminants file:
1. Sequences must be on a single line (ie can't contain line breaks)
2. No blank lines in file.
3. Sequence header is space delimited list of options. Option names and values should be separated by a colon. 
   Ex "> name:oligo size:10 windows:5"
4. Options:
  name: this option is required for every sequence
  method: there are three trimming methods (0 = full contaminant, 1 = mapped contaminant, 2 = identity based). (default 2)
  size: size of k-mer (default 7)
  windows: how many k-mers to seek. can not be larger than (contaminant length - k-mer). (default 6)
  percentIdentity: percent identity threshold for trimming method 2 (0.0 < percentIdentity <= 1.0). (default 0.9)
  totalIdentity: total identity threshold for trimming method 2. If this is less than the k-mer size then it will have no impact on trimming. (default 16)
"""

#------------------------------------------------------------------------------------------
# INITIALIZATIONS
#------------------------------------------------------------------------------------------

import sys, os, argparse

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

VERSION = 0.4

# indecies for the read set
HEADER = 'header'
SEQUENCE = 'seq'
QUALS = 'quals'
LENGTH = 'length'
TRIMMED = 'trimmed'

DEFAULT_END = 3
DEFAULT_MAX = 0
DEFAULT_METHOD = 2
DEFAULT_PERCENT_IDENTITY = .9
DEFAULT_SIZE_KMER = 7
DEFAULT_TOTAL_IDENTITY = 16
DEFAULT_WINDOWS_NUMBER = 6

# used by removePolyAT()
POLY_A = 'AAAAAAAAAA'
POLY_T = 'TTTTTTTTTT'

argParser = argparse.ArgumentParser(version=str(VERSION), 
                                    description='Trim NGS reads.',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    epilog='' +
'Trimming will happen in the following order depending on which options are selected:\n' + 
'\t1. Cut 3\' end of read.\n' + 
'\t2. Cut 5\' end of read.\n' + 
'\t3. Remove N\'s from both ends.\n' +
'\t4. Process contaminants file, removing contaminants based on their order in the contaminants file.\n' +
'\t5a. Single-end: discard read if shorter than the minimum length.\n' + 
'\t5b. Paired-end: if only one of the paired reads is shorter than the minimum length, then replace that read\'s sequence with N\'s and replace that read\'s quality scores with #. If both paired reads are shorter than the minimum length, then discard the read pair.\n' + 
'\t6. Pad paired reads with N\'s so that they are both the same length. For every N that is added, also add a # to the quality score.\n' + 
'\t7. Add "L:X" to the read header with X being the new length of the sequence (including any N\'s that were added to the sequence).\n\n' + 

'Contaminants file:\n' + 
'\t1. Sequences must be on a single line (ie can\'t contain line breaks).\n' +
'\t2. No blank lines in file.\n' +
'\t3. Sequence header is space delimited list of options. Option names and values should be separated by a colon. \n' +
'\t\t Example header "> name:oligo end:3 size:10 windows:5"\n' +
'\t4. Options:\n' +
'\t\t * name: this option is required for every sequence\n' +
'\t\t * method: there are three trimming methods (0 = full contaminant, 1 = mapped contaminant, 2 = identity based). (default 2)\n' +
'\t\t\t 0. Full contaminant trimming means that when a k-mer is mapped then it is expected that the entire contaminant mapped and the read is trimmed accordingly. For example lets assume we have a k-mer that is located 4 bases from the 5\' end of a contaminant. If that k-mer maps then we would shift where we trim the read by 4 bases in the direction of the 5\' end of the read. We would then remove all bases from that position to the 3\' end of the read, regardless of the additional bases mapped to the contaminant.\n' +
'\t\t\t 1. Mapped contaminant trimming means that when a k-mer is mapped then we extend the mapping and trimmed accord to the mapping. For example lets assume we have a k-mer that is located 4 bases from the 5\' end of a contaminant. If that k-mer maps then we would extend the mapped region one base at a time, in the 5\' direction until we found a base that didn\'t map. We would then trim from that postion to the 3\' end of the read.\n' +
'\t\t\t 2. If a k-mer maps to the read then the location of the mapping is used to anchor the contaminant to the read. The percent and total identity between the contaminant and the read is computed. If both the percent and total identity are above a user-defined threshold then the read is trimmed from the beginning of the contaminant to the 3\' end of the read. If not then the read is not trimmed.\n' +
'\t\t * size: size of k-mer (default 7)\n' +
'\t\t * windows: how many k-mers to seek. can not be larger than (contaminant length - k-mer). (default 6)\n' +
'\t\t * percentIdentity: percent identity threshold for trimming method 2 (0.0 < percentIdentity <= 1.0). (default 0.9)\n' +
'\t\t * totalIdentity: total identity threshold for trimming method 2. If this is less than the k-mer size then it will have no impact on trimming. (default 16)\n'
                                    )
argParser.add_argument( '-p', '--padPaired', dest='padPaired', action='store_true', default=False,
                        help='Pad paired reads so that they are the same length after trimming all trimming has occured. N\'s will be added to the 3\' end with \'#\' added to the quality score for each N that is added. This will not do anything for single-end reads. (default: no)' )
argParser.add_argument( '-m', '--minLen', dest='min_len', action='store', default=0, type=int,
                        help='Minimum size of trimmed read. If trimmed beyond minLen, then read is discarded. If read is paired then read is replaced with N\'s, unless both reads in pair are smaller than minLen in which case the pair is discarded. (default: no minimum length)' )
argParser.add_argument( '-c3', '--cut3', dest='num_cut_3', action='store', default=0, type=int, 
                        help='number of bases to remove from 3\' end of read. Truncating reads does not count toward trimming totals. This happens prior to the removing of N\'s and hence prior to contaminant trimming. (default: 0)' )
argParser.add_argument( '-c5', '--cut5', dest='num_cut_5', action='store', default=0, type=int, 
                        help='number of bases to remove from 5\' end of read. Truncating reads does not count toward trimming totals. This happens prior to the removing of N\'s and hence prior to contaminant trimming. (default: 0)' )
argParser.add_argument( '-rN', '--removeNs', dest='removeN', action='store_true', default=False,
                        help='remove N\'s from both ends of the read. This trimming happens before contaminant trimming. (default: no)' )
argParser.add_argument( '-rAT', '--removePolyAT', dest='remove_AT', action='store', default=-1, type=int,
                        help='length of 3\' poly-A and 5\' poly-T to remove from the respective ends of the read. If all poly A/T is to be removed then the value should be equal to or greater than the length of the read. A minimum of ten A\'s or ten T\'s must exist in order for this trimming to happen, regardless of the trimming length; that is, poly-A and poly-T fragments are defined as being at least 10 nt in length. A sequences of A\'s or T\'s are ignored. This trimming happens after contaminant trimming. (default: no trimming)' )
argParser.add_argument( '-c', dest='contaminants_fa', action='store', default=None, 
                        help='fasta-like file containing list of contaminants to trim from the 3\' end of the read' )
argParser.add_argument( '-f', dest='forward_fq', action='store', required=True,
                        help='fastq file with reads to be trimmed' )
argParser.add_argument( '-r', dest='reverse_fq', 
                        help='second read with paired-end reads. This file is not present for single-end reads.' )
argParser.add_argument( '-o', dest='output_prefix', action='store', required=True,
                        help='prefix for output file(s). A \'_1\' will be appended to the forward reads output file and if paired reads then a \'_2\' will be appended to the reverse reads file. Output files similarly named and with the suffix \'.disc.txt\' will be created to store the fastq headers for reads shorter than the minimum length threshold after trimming.' )

clArgs = argParser.parse_args()
if DEBUG: print clArgs

# flag if paired-end reads
PAIRED = False
if clArgs.reverse_fq:
    PAIRED = True

# track trimming stats
nBothTrimmed = 0 # total number of read pairs in which both reads were trimmed (equal to nForTrimmed if single-end)
nForTrimmed = 0 # num forward reads trimmed
nRevTrimmed = 0 # num reverse reads trimmed
nBothDiscarded = 0 # num read pairs discarded (equal to nForDiscarded if single-end)
nForDiscarded = 0 # num of forward reads replaced with N's
nRevDiscarded = 0 # num of reverse reads replaced with N's
nTotalReadPairs = 0 # total number of read pairs processed
nForContaminantsTrim = {} # number of each contaminant trimmed from forward read
nRevContaminantsTrim = {} # number of each contaminant trimmed from reverse read

# print error and quit
def quitOnError(msg):
    print 'ERROR:', msg
    sys.exit(0)

# generate output file with every trim event. Forward and reverse
# reads will be combined in this output file.
if DEBUG:
    try: 
        debugFile = open("debug.tsv", 'w')
        debugFile.write('Orig Seq\tLen\tTrimmed Seq\tLen\tTrim Type\tMisc\n')
    except: 
        quitOnError('Unable to open output file debug.tsv')

#------------------------------------------------------------------------------------------
# FUNCTIONS FOR HANDLING READING AND WRITING OF READS
#------------------------------------------------------------------------------------------

def nextRead(inFile):
    """
    Returns a set consisting of each line from the read, a length
    indicator and a trimmed flag. Returns empty set if no more reads.
    """
    line = inFile.readline().rstrip()
    if len(line) == 0:
        return {}

    read = {}
    try:
        # load read values from fastq file
        read[HEADER] = line
        # make sure sequence is upper case for comparisons
        read[SEQUENCE] = inFile.readline().rstrip().upper()
        # this is the + which we can ignore for now

        # XXX CONFIRM THAT THE + IN THE FASTQ FILE IS A CONSTANT.
        inFile.readline()

        read[QUALS] = inFile.readline().rstrip()
        
        # store additional read information in read set
        read[LENGTH] = len(read[SEQUENCE])
        read[TRIMMED] = False

        return read
    except:
        msg = 'Problem loading read:', line
        quitOnError(msg)

def writeRead(read, outFile):
    """
    Writes a read to the read file. Read is expected to be a 4-item
    list as returned by nextRead().
    """
    outFile.write(read[HEADER] + '\n') # header
    outFile.write(read[SEQUENCE] + '\n') # sequence
    outFile.write('+\n') # +
    outFile.write(read[QUALS] + '\n') # quals

def debugTrimOutput(origSeq, origLen, trimSeq, trimLen, trimType, misc):
    """
    Output information about the trim event.
    """
    debugFile.write(origSeq + '\t' + str(origLen) + '\t' + trimSeq + '\t' + str(trimLen) + '\t' + trimType + '\t' + misc + '\n')
    
#------------------------------------------------------------------------------------------
# TRIMMING FUNCTIONS
#------------------------------------------------------------------------------------------

# truncate 3' end of reads. Truncated reads don't count as having been
# trimmed as this impacts all reads.
def cut3(read, nCut):
    """
    Cuts nCut bases from 3' end of read.
    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    length = read[LENGTH]

    # trim sequence
    seq = seq[nCut:]

    # need to trim quals in same way we trimmed sequence
    quals = quals[nCut:]

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    return read

# truncate 5' end of reads. Truncated reads don't count as having been
# trimmed as this impacts all reads.
def cut5(read, nCut):
    """
    Cuts nCut bases from 5' end of read.
    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    length = read[LENGTH]

    # trim sequence
    seq = seq[:(length - nCut)]

    # need to trim quals in same way we trimmed sequence
    quals = quals[:(length - nCut)]

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    return read

# remove N's on either end of the read
def removeNs(read):
    """
    Removes any N's from the 3' and 5' ends of a sequence.
    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    length = read[LENGTH]

    # local flag for trimming having occurred. This is unnecessary
    # since removeNs() is first trimming to happen. However this
    # prevents possible problem in future if we add another trimming
    # option prior to removeNs()
    wasTrimmed = False 

    # trim N from beginning of read (5' end)
    if seq.startswith('N'):
        # trim sequence
        seq = seq.lstrip('N')

        # need to trim quals in same way we trimmed sequence
        quals = quals[len(quals) - len(seq):]

        # flag read as having been trimmed
        wasTrimmed = True

    # trim N from end of read (3' end)
    if seq.endswith('N'):
        # trim sequence
        seq = seq.rstrip('N')

        # need to trim quals in same way we trimmed sequence
        quals = quals[:len(seq)]

        # flag read as having been trimmed
        wasTrimmed = True

    if wasTrimmed:
        if DEBUG: debugTrimOutput(read[SEQUENCE], length, seq, len(seq), 'removeNs', '')

        read[SEQUENCE] = seq
        read[QUALS] = quals
        read[LENGTH] = len(seq)
        read[TRIMMED] = True
    return read

# trim contaminant removing entire contaminant based on single k-mer
# mapping, from the 5' end of the contaminant to the 3' end of the
# read.
def fullContaminantTrimming(read, contaminant, isForRead):
    """
    Trim read based on contaminant. Contaminant is expected to be a
    tuple containing a set of trimming options and the contaminant
    sequence. For trimming first map a k-mer then trim the read from
    the expected start point of the contaminant based on the position
    of the k-mer within the contaminant to the end of the read. Note
    that this function is an agressive trimming scheme and may remove
    portions of the read that don't map to the contaminant.
    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    length = read[LENGTH]

    # list of contaminant k-mers to look for in read
    kmers = contaminant[2]
    # load contaminant's options
    options = contaminant[0]
    # contaminant name
    name = options['name']

    # we attempt to perfectly align each k-mer against the
    # sequence. If any alignment succeeds, the position of the 5' end
    # of the contaminant is used for trimming. For example if the
    # mapped k-mer begins at the 5th base in sequence A then sequence
    # A is assumed to extend another 4 bases toward the 5' end of
    # sequence B and hence that position in sequence B is used for
    # trimming. If none of the k-mer's map then no trimming occurs.
    pos = -1
    for kmer, offset in kmers:
        # look for k-mer in read sequence. Index is the 5' position in
        # the read where the k-mer mapped. If the k-mer maps multiple
        # times in the read, then the mapping closest to the 3' end is
        # used.
        index = seq.rfind(kmer)
        if index > -1: 
            # the read sequence may not contain the entire contaminant
            # so change the order of the k-mer search based on the end
            # we are searching. Search the k-mers that are located
            # near the 5' end of the contaminant before searching for
            # the k-mers located near the 3' end of the contaminant.
            pos = index - offset

            # break out of loop since we found a match
            break

    # if k-mers not found, then no trimming
    if pos == -1: return read

    seq = seq[:pos]
    quals = quals[:pos]

    if DEBUG: debugTrimOutput(read[SEQUENCE], length, seq, len(seq), 'fullContaminantTrimming', name)

    if isForRead: nForContaminantsTrim[name] += 1
    else: nRevContaminantsTrim[name] += 1

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    read[TRIMMED] = True
    return read

# trim contaminant removing the portion of the contaminant that maps
# from the k-mer region to the end of the sequence.
def mappedContaminantTrimming(read, contaminant, isForRead):
    """Trim read based on contaminant. Contaminant is expected to be a
    tuple containing a set of trimming options and the contaminant
    sequence. For trimming first map a k-mer then extend the mapping
    toward the 5' end of the read. After the mapping is fully
    extended, trim the read from the contaminant mapping region to the
    3' end of the read. Note that this function is a more conservative
    trimming scheme and may miss contaminate regions that are smaller
    than the size of a k-mer.

    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    length = read[LENGTH]

    # contaminant sequence
    cSeq = contaminant[1]
    cLength = len(cSeq)
    # list of contaminant k-mers to look for in read
    kmers = contaminant[2]
    # load contaminant's options
    options = contaminant[0]
    # contaminant name
    name = options['name']

    # we attempt to perfectly align each k-mer against the
    # sequence. Once we find a match, then we extend the mapping
    # toward the 5' end of the read, only allowing for exact
    # matches. We then trim from the mapped region of the read to the
    # 3' end of the read. For example if the mapped k-mer begins at
    # the 5th base in seuence A then we extend the k-mer mapped toward
    # the 5' end of the sequences. The position returned is the last
    # mapped base on the 5' end of the mapping. If none of the k-mer's
    # map then no trimming occurs.
    pos = -1
    for kmer, offset in kmers:
        # when we set pos we break out of the inner loop. Here we
        # break out of the outer loop.
        if pos > -1: break

        # look for k-mer in read sequence. Index is the 5' position in
        # the read where the k-mer mapped. We're expecting the k-mer
        # to map to the 3' end of the read so we search the k-mers
        # that are located near the 5' end of the contaminant before
        # searching for the k-mers located near the 3' end of the
        # contaminant.  If the k-mer maps multiple times in the read,
        # then the mapping closest to the 3' end is used.
        index = seq.rfind(kmer)
        if index > -1: 
            # found mapping so now extend mapping toward 5' end

            pos = 0
            while offset > 0 and index > 0:
                # offset is the position of the 5' end of the k-mer
                # in the contaminant. So if equal to 0 then we are
                # already at the 5' end of the contaminant. index
                # is our position in the read sequence. If that's
                # 0 then we are at the 5' end of the read.
                if seq[index-1] != cSeq[offset-1]:
                    pos = index
                    break
                else:
                    index -= 1
                    offset -= 1

    # if k-mers not found, then no trimming
    if pos == -1: return read

    seq = seq[:pos]
    quals = quals[:pos]

    if DEBUG: debugTrimOutput(read[SEQUENCE], length, seq, len(seq), 'mappedContaminantTrimming', name)

    if isForRead: nForContaminantsTrim[name] += 1
    else: nRevContaminantsTrim[name] += 1

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    read[TRIMMED] = True
    return read

# trim contaminant removing entire contaminant based on single k-mer
# mapping and whether or not the percent identity between the
# contaminant and the read is above the specified thresholds
def identityTrimming(read, contaminant, isForRead):
    """
    Trim read based on contaminant. Contaminant is expected to be a
    tuple containing a set of trimming options and the contaminant
    sequence. For trimming first map a k-mer to align the contaminant
    to the read. Then compute the percent identity between the
    contaminant and the read. If the percent identity is above the
    specified threshold then trim read from the position of the 5' end
    of the contaminant to the 3' end of the read. If the threshold is
    not reached then do not trim the read. With a high threshold this
    function will bias away from trimming. As the threshold is
    decreased (toward 0) then the likelyhood of trimming increases.
    Note that the identity is only computed for the portion of the
    read that overlaps with the contaminant. So if a short contaminant
    perfectly maps in the middle of a long read then it would have a
    high percent identity, causing the read to be trimmed.
    Contaminant example: [{'name': 'indexAdapter', 'windows': 6,
    'percentIdentity': 0.9, 'method': 2, 'pecentIdentity': 0.9,
    'totalIdentity': 16, 'size': 10},
    'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', [('AGATCGGAAG', 0),
    ('CGGAAGAGCA', 4), ('AGAGCACACG', 8), ('CACACGTCTG', 12),
    ('CGTCTGAACT', 16), ('TGAACTCCAG', 20), ('CTCCAGTCAC', 24)]]
    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    length = read[LENGTH]

    # contaminant sequence
    cSeq = contaminant[1]
    cLength = len(cSeq)
    # list of contaminant k-mers to look for in read
    kmers = contaminant[2]
    # load contaminant's options
    options = contaminant[0]
    # contaminant name
    name = options['name']
    # percent identity threshold
    percIdentityThreshold = options['percentIdentity']
    # total identity threshold
    totIdentityThreshold = options['totalIdentity']

    # we attempt to perfectly align each k-mer against the
    # sequence. If any alignment succeeds, then we compute the percent
    # identity between the contaminant and read. If above our
    # threshold then trim the entire contaminant, otherwise do not
    # trim anything.
    percIdentity = 0
    totIdentity = 0
    for kmer, offset in kmers:
        # look for k-mer in read sequence. Index is the position from
        # the 5' end of the read where the k-mer mapped to the read.
        # If the k-mer maps multiple times in the read, then the
        # mapping closest to the 3' end is used.
        index = seq.rfind(kmer)
        if index > -1: 
            # Pos is the position where the 5' end of contaminant
            # overlaps the read, based in k-mer mapping. Offset is
            # the position of the 5' end of the k-mer within the
            # contaminant.
            pos = index - offset

            # cPos is the position in the contaminant
            cPos = 0

            # if the contaminant hangs off the 5' end of the read then
            # we need to adjust the search to begin from the 5' end of
            # the read
            if pos < 0:
                # shift where we begin in the contaminant based on how
                # many bases are hanging off the 5' end of the read
                cPos = pos * -1
                # start at 5' end of read
                pos = 0

            # rPos is the position in the read
            rPos = pos

            totIdentity = 0
            count = 0
            while cPos < cLength and rPos < length:
                # count the number of matching bases
                if seq[rPos] == cSeq[cPos]: totIdentity += 1
                
                # increment positions so we check next base
                cPos += 1
                rPos += 1
                
                # count the number of bases checked
                count += 1

            # we've now checked all overlapping bases, so we can
            # compute the percent identity.
            percIdentity = float(totIdentity) / float(count)

            # break out of loop since we found a match
            break

    # if k-mers not found, then no trimming
    if index == -1: return read

    # if the k-mers were not found or the contaminant didn't map
    # sufficiently to the read then we don't trim.
    if (percIdentity < percIdentityThreshold) or (totIdentity < totIdentityThreshold): return read

    seq = seq[:pos]
    quals = quals[:pos]

    if DEBUG: 
        misc = name + "\t" + str(percIdentity) + "\t" + str(totIdentity)
        debugTrimOutput(read[SEQUENCE], length, seq, len(seq), 'identityTrimming', misc)

    if isForRead: nForContaminantsTrim[name] += 1
    else: nRevContaminantsTrim[name] += 1

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    read[TRIMMED] = True
    return read

# remove poly A from 3' end and poly T from 5' end
def removePolyAT(read, trimLength):
    """
    Removes trimLength number of A's from the 3' end and T's from the
    5' end of a sequence. First look to see if there are at least 10
    A's or 10 T's. If so then trim all A's or all T's.
    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    length = read[LENGTH]

    wasTrimmed = False # local flag for trimming

    # trim poly-T from 5' of read
    if seq.startswith(POLY_T):
        # only remove up to trimLength number of bases from poly-T.
        # Get copy of trimLength number of bases in sequence.
        n_mer = seq[:trimLength]
        # remove poly-T from 5' end of n-mer
        n_mer = n_mer.lstrip('T')
        # append remaining bases, if any, back onto 5' end of read
        seq = n_mer + seq[trimLength:]

        # need to trim quals to same degree that we trimmed sequence
        quals = quals[len(quals) - len(seq):]

        # flag read as having been trimmed
        wasTrimmed = True

    # trim poly-A from 3' end of read
    if seq.endswith(POLY_A):
        # only remove up to trimLength number of bases from poly-A.
        # Get copy of last trimLength number of bases in sequence.
        idx = len(seq) - trimLength
        n_mer = seq[idx:]
        # remove poly-A from 3' end of n-mer
        n_mer = n_mer.rstrip('A')
        # append remaining bases back onto 3' end of read
        seq = seq[:idx] + n_mer

        # need to trim quals in same way we trimmed sequence
        quals = quals[:len(seq)]

        # flag read as having been trimmed
        wasTrimmed = True

    if wasTrimmed:
        if DEBUG: debugTrimOutput(read[SEQUENCE], length, seq, len(seq), 'removePolyAT', '')

        read[SEQUENCE] = seq
        read[QUALS] = quals
        read[LENGTH] = len(seq)
        read[TRIMMED] = True
    return read

#------------------------------------------------------------------------------------------
# LOAD CONTAMINANTS
#------------------------------------------------------------------------------------------

def computeKMers(seq, kmerSize, numWindows):
    """
    This will return a tuple of k-mers that includes the k-mer
    sequence. The tuple will also include either 5' position of the
    k-mer within the original sequence. The k-mers are ordered from
    the 5' end of the contaminant to the 3' end. We don't expect the
    read sequence to contain the entire contaminant so we start at the
    5' end when we are searching. Expecting the k-mer to map to the 3'
    end of the read then we search the k-mers that are located near
    the 5' end of the contaminant before searching for the k-mers
    located near the 3' end of the contaminant.
    """
    seqLength = len(seq)

    # the last window is bounded by the length of the sequence and the size of a k-mer
    lastWindow = seqLength - kmerSize

    # divide up the length of the sequence allowing for the first
    # k-mer to be at the beginning of the sequence and the last k-mer
    # to be approximate at the end of the sequence. The ordering of
    # the k-mers is dependent on which end of the read is expected to
    # contain the contaminant. Depending on the length of the
    # contaminant, the k-mer size and the number of windows, the last
    # few bases on the tail end of the contaminant might not be
    # included in any k-mer. For a contaminant expected to be found on
    # the 3' end of a read this means the last few bases on the 3' end
    # of the contaminant might not be in any k-mer. For a contaminant
    # expected to be found on the 5' end of a read then the last few
    # bases on the 5' end of the contaminant might not be in a k-mer.
    if numWindows > 1:
        windowStep = lastWindow / (numWindows - 1)
    else:
        # make sure we only get one k-mer
        windowStep = seqLength

    kmerList = []
    # slice up the sequence based on the window step size and the k-mer size
    for index in range(0, lastWindow + 1, windowStep):
        # store k-mer and the k-mer position offset for a 3'
        # mapping. Offset is the 5' position of the k-mer.
        kmerList.append((seq[index:index+kmerSize], index))

    return kmerList

if clArgs.contaminants_fa:
    try: contFile = open(clArgs.contaminants_fa, 'r')
    except IOError: print 'ERROR: No such file', clArgs.contaminants_fa

    # need to retain order so using tuple instead of set
    contaminantList = []
    count = 0
    for line in contFile:
        if line.startswith('>'):
            header = line[1:].strip()
            # need to process header arguments
            args = header.split()

            # convert options string into set. Initialize values to
            # defaults. The initialized set doesn't include the name
            # as name is required.
            options = { 'method':DEFAULT_METHOD, 
                        'size':DEFAULT_SIZE_KMER, 
                        'windows':DEFAULT_WINDOWS_NUMBER,
                        'percentIdentity':DEFAULT_PERCENT_IDENTITY, 
                        'totalIdentity':DEFAULT_TOTAL_IDENTITY }
            for option in args:
                key, value = option.split(':')

                if key == 'method':
                    value = int(value)
                    if value != 0 and value != 1 and value != 2:
                        msg = 'Invalid "method" option (%d) for contaminant. Must be 0, 1, or 2.' % value
                        quitOnError(msg)
                elif key == 'size':
                    value = int(value)
                elif key == 'windows':
                    value = int(value)
                elif key == 'totalIdentity':
                    value = int(value)
                elif key == 'percentIdentity':
                    value = float(value)
                    if value <= 0.0 or value > 1.0:
                        msg = 'Invalid "percentIdentity" option (%f) for contaminant. Percent identity must be greater than 0 and less than or equal to 1.' % value
                        quitOnError(msg)
                # at this point we've accounted for all options but
                # the name. If we don't have a name then the option
                # doesn't exist.
                elif key != 'name':
                    msg = 'Invalid contamination option "%s".' % option
                    quitOnError(msg)

                # save value in options
                options[key] = value

            # make sure we loaded a name.
            if 'name' not in options:
                quitOnError('Contaminant does not have a name.')

            # use name to initialize contaminant trimming counts
            nForContaminantsTrim[options['name']] = 0
            if PAIRED: nRevContaminantsTrim[options['name']] = 0
        else:
            seq = line.strip().upper()

            # compute set of k-mers
            kmerList = computeKMers(seq, options['size'], options['windows'])

            if len(seq) < options['size']:
                msg = 'The k-mer size (%d) must be smaller than the length of the sequence (%s).' % (options['size'], seq)
                quitOnError(msg)

            # save contaminant and related values in tuple
            contaminantList.append([options, seq, kmerList])

            count += 1

    contFile.close()
    
    print 'Loaded contaminants: ' + str(count)
    if DEBUG: 
        for c in contaminantList: print c

#------------------------------------------------------------------------------------------
# OPEN READ INPUT AND OUTPUT FILES
#------------------------------------------------------------------------------------------

# open input file
forReadIn = ''
try: 
    forReadIn = open(clArgs.forward_fq, 'r')
except: 
    msg = 'Unable to load reads file ' + clArgs.forward_fq
    quitOnError(msg)

# open file that will store non-discarded reads
forReadOut = ''
try: 
    forReadOut = open(clArgs.output_prefix + "_1.fq", 'w')
except: 
    msg = 'Unable to open output file ' + clArgs.output_prefix + "_1.fq"
    quitOnError(msg)

# open file that will store discarded reads
forReadDiscardedOut = ''
try: 
    forReadDiscardedOut = open(clArgs.output_prefix + "_1.disc.txt", 'w')
except: 
    msg = 'Unable to open output file ' + clArgs.output_prefix + "_1.disc.txt"
    quitOnError(msg)

if PAIRED:
    revReadIn = ''
    try: 
        revReadIn = open(clArgs.reverse_fq, 'r')
    except: 
        msg = 'Unable to load reads file ' + clArgs.reverse_fq
        quitOnError(msg)

    revReadOut = ''
    try: 
        revReadOut = open(clArgs.output_prefix + "_2.fq", 'w')
    except: 
        msg = 'Unable to open output file ' + clArgs.output_prefix + "_2.fq"
        quitOnError(msg)

    revReadDiscardedOut = ''
    try: 
        revReadDiscardedOut = open(clArgs.output_prefix + "_2.disc.txt", 'w')
    except: 
        msg = 'Unable to open output file ' + clArgs.output_prefix + "_2.disc.txt"
        quitOnError(msg)

#------------------------------------------------------------------------------------------
# TRIM READS.  LOAD ONE READ FROM BOTH FQ FILES AT SAME TIME. PROCESS
# READ PAIR TOGETHER AND THEN EITHER WRITE TO OUTPUT FILE OR DISCARD.
# ------------------------------------------------------------------------------------------

while 1:
    forRead = nextRead(forReadIn)
    if PAIRED: revRead = nextRead(revReadIn)

    if forRead == {}:
        forReadIn.close()
        forReadOut.close()
        if PAIRED: 
            revReadIn.close()
            revReadOut.close()
        break
            
    nTotalReadPairs += 1

    # use flags since we only want to count once every time a read is trimmed
    forReadTrimmed = False
    revReadTrimmed = False

    #--------------------------------------------------------------------------------------
    # cut 3' end
    if clArgs.num_cut_3 > 0:
        forRead = cut3(forRead, clArgs.num_cut_3)
        if PAIRED: revRead = cut3(revRead, clArgs.num_cut_3)

    #--------------------------------------------------------------------------------------
    # cut 5' end
    if clArgs.num_cut_5 > 0:
        forRead = cut5(forRead, clArgs.num_cut_5)
        if PAIRED: revRead = cut5(revRead, clArgs.num_cut_5)

    #--------------------------------------------------------------------------------------
    # remove N's
    if clArgs.removeN:
        forRead = removeNs(forRead)
        if PAIRED: revRead = removeNs(revRead)

    #--------------------------------------------------------------------------------------
    # trim contaminants
    if clArgs.contaminants_fa:
        for contaminant in contaminantList:
            method = contaminant[0]['method']
            if method == 0:
                forRead = fullContaminantTrimming(forRead, contaminant, True)
                if PAIRED: revRead = fullContaminantTrimming(revRead, contaminant, False)
            elif method == 1:
                forRead = mappedContaminantTrimming(forRead, contaminant, True)
                if PAIRED: revRead = mappedContaminantTrimming(revRead, contaminant, False)
            else:
                forRead = identityTrimming(forRead, contaminant, True)
                if PAIRED: revRead = identityTrimming(revRead, contaminant, False)

    #--------------------------------------------------------------------------------------
    # remove poly A/T
    if clArgs.remove_AT > 0:
        forRead = removePolyAT(forRead, clArgs.remove_AT)
        if PAIRED: revRead = removePolyAT(revRead, clArgs.remove_AT)

    #--------------------------------------------------------------------------------------
    # compute trimming stats
    trimFor = False
    trimRev = False
    if forRead[TRIMMED]: 
        nForTrimmed += 1
        trimFor = True
    if PAIRED: 
        if revRead[TRIMMED]: 
            nRevTrimmed += 1
            trimRev = True
        if trimFor and trimRev:
            # count trimmed pair
            nBothTrimmed += 1
    elif trimFor:
        # not paired, so count single read as pair
        nBothTrimmed += 1

    #--------------------------------------------------------------------------------------
    # discard reads that are too short
    if clArgs.min_len > 0:
        discardFor = False
        discardRev = False
        if forRead[LENGTH] < clArgs.min_len:
            nForDiscarded += 1
            discardFor = True
        if PAIRED:
            if revRead[LENGTH] < clArgs.min_len:
                nRevDiscarded += 1
                discardRev = True
            if discardFor and discardRev:  
                # count pair
                nBothDiscarded += 1

                # both reads in read pair are discarded so don't write
                # them to output file. Instead write the read headers
                # to the discard file. We only store the headers in
                # the file since, in some cases, the entire sequence
                # will have been trimmed and hence this would be an
                # invalid fastq file.
                forReadDiscardedOut.write(forRead[HEADER] + '\n') # header
                revReadDiscardedOut.write(revRead[HEADER] + '\n') # header
                continue
        elif discardFor:
            # not paired so count single read as pair
            nBothDiscarded += 1

            # single-end read being discarded, so don't write to
            # output file. Instead write the read headers to the
            # discard file.
            forReadDiscardedOut.write(forRead[HEADER] + '\n') # header
            continue

        # if we got this far then paired reads with only one being too
        # small, so replace that one read with N's
        if discardFor:
            forRead[SEQUENCE] = 'N' * revRead[LENGTH]
            forRead[QUALS] = '#' * revRead[LENGTH]
            forRead[LENGTH] = revRead[LENGTH]
        elif discardRev:
            revRead[SEQUENCE] = 'N' * forRead[LENGTH]
            revRead[QUALS] = '#' * forRead[LENGTH]
            revRead[LENGTH] = forRead[LENGTH]

    #--------------------------------------------------------------------------------------
    # pad paired reads, as needed, so that they are both the same length
    if clArgs.padPaired and PAIRED:
        if forRead[LENGTH] < revRead[LENGTH]:
            forRead[SEQUENCE] = forRead[SEQUENCE] + 'N' * (revRead[LENGTH] - forRead[LENGTH])
            forRead[QUALS] = forRead[QUALS] + '#' * (revRead[LENGTH] - forRead[LENGTH])
            forRead[LENGTH] = revRead[LENGTH]
        elif revRead[LENGTH] < forRead[LENGTH]:
            revRead[SEQUENCE] = revRead[SEQUENCE] + 'N' * (forRead[LENGTH] - revRead[LENGTH])
            revRead[QUALS] = revRead[QUALS] + '#' * (forRead[LENGTH] - revRead[LENGTH])
            revRead[LENGTH] = forRead[LENGTH]

    #--------------------------------------------------------------------------------------
    # add sequence length to read header
    forRead[HEADER] += ' L:%d' % forRead[LENGTH]
    if PAIRED: revRead[HEADER] += ' L:%d' % revRead[LENGTH]

    #--------------------------------------------------------------------------------------
    # write read(s) to output file if not over-trimmed
    writeRead(forRead, forReadOut)
    if PAIRED: writeRead(revRead, revReadOut)

#------------------------------------------------------------------------------------------
# OUTPUT TRIMMING STATS
#------------------------------------------------------------------------------------------

print 'Read pairs processed:', nTotalReadPairs
print 'Both forward and reverse reads trimmed:', nBothTrimmed
print '\tForward reads trimmed:', nForTrimmed
print '\tReverse reads trimmed:', nRevTrimmed
print 'Both forward and reverse reads discarded:', nBothDiscarded
print '\tForward reads discarded:', nForDiscarded
print '\tReverse reads discarded:', nRevDiscarded
print 
for name, count in nForContaminantsTrim.iteritems() :
    print 'Contaminant: ', name
    print '\tForward reads trimmed', count
    nRev = 0
    if PAIRED: nRev = nRevContaminantsTrim[name]
    print '\tReverse reads trimmed', nRev

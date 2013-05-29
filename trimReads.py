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

usage: trimReads.py [-h] [-v] [-p] [-m MINLEN] [-c3 NUM_CUT_3] [-c5 NUM_CUT_5] [-rN] [-c CONTAMINANTS_FA] -f
                    FORWARD_FQ [-r REVERSE_FQ] -o OUTPUT_PREFIX

Contaminants file:
1. Sequences must be on a single line (ie can't contain line breaks)
2. No blank lines in file.
3. Sequence header is space delimited list of options. Option names and values should be separated by a colon. 
   Ex "> name:oligo end:3 size:10 windows:5"
4. Options:
  end: expect contaminant on the 3' or 5' end. (default 3)
  max: maximum number of bases to trim. useful with poly A/T trimming. (default all matches)
  method: there are three trimming methods (0 = full contaminant, 1 = mapped contaminant, 2 = identity based). (default 2)
  name: this option is required for every sequence
  percentIdentity: percent identity threshold for trimming method 2 (0.0 < percentIdentity <= 1.0). (default 0.9)
  size: size of k-mer (default 7)
  totalIdentity: total identity threshold for trimming method 2. If this is less than the k-mer size then it will have no impact on trimming. (default 16)
  windows: how many k-mers to seek. can not be larger than (contaminant length - k-mer). (default 6)

**********************************************************************
TODO
- add option to output discarded reads to separate output file
- use contaminant name to print contaminant specific stats
- need to add minimum mapping length to method 2. we need to check the total identity along with the percent identity.
**********************************************************************

"""

#------------------------------------------------------------------------------------------
# INITIALIZATIONS
#------------------------------------------------------------------------------------------

import sys, os, argparse

DEBUG = True
if DEBUG: print 'DEBUG MODE: ON'

VERSION = 0.1

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
'\t\t * end: expect contaminant on the 3\' or 5\' end (values: 3 or 5). (default 3)\n' +
'\t\t * max: maximum number of bases to trim. useful with poly A/T trimming. (default all matches)\n' +
'\t\t * method: there are three trimming methods (0 = full contaminant, 1 = mapped contaminant, 2 = identity based). (default 2)\n' +
'\t\t\t 0. Full contaminant trimming means that when a k-mer is mapped then it is expected that the entire contaminant mapped and the read is trimmed accordingly. For example lets assume we have a k-mer that is located 4 bases from the 5\' end of a contaminant. Furthermore lets assume the contaminant is expected to be located on the 3\' end if a read. If that k-mer maps then we would shift where we trim the read by 4 bases in the direction of the 5\' end of the read. We would then remove all bases from that position to the 3\' end of the read, regardless if the additional bases mapped to the contaminant.\n' +
'\t\t\t 1. Mapped contaminant trimming means that when a k-mer is mapped then we extend the mapping and  trimmed accord to the mapping. For example lets assume we have a k-mer that is located 4 bases from the 5\' end of a contaminant. Furthermore lets assume the contaminant is expected to be located on the 3\' end if a read. If that k-mer maps then we would extend the mapped region one base at a time, in the 5\' direction until we found a base that didn\'t map. We would then trim from that postion to the 3\' end of the read.\n' +
'\t\t\t 2. If a k-mer maps to the read then the location of the mapping is used to anchor the contaminant to the read. The percent and total identity between the contaminant and the read is computed. If both the percent and total identity are above a user-defined threshold then the read is trimmed from the beginning of the contaminant to the end of the read. If not then the read is not trimmed.\n' +
'\t\t * name: this option is required for every sequence\n' +
'\t\t * percentIdentity: percent identity threshold for trimming method 2 (0.0 < percentIdentity <= 1.0). (default 0.9)\n' +
'\t\t * size: size of k-mer (default 7)\n' +
'\t\t * totalIdentity: total identity threshold for trimming method 2. If this is less than the k-mer size then it will have no impact on trimming. (default 16)\n' +
'\t\t * windows: how many k-mers to seek. can not be larger than (contaminant length - k-mer). (default 6)\n'
                                    )
argParser.add_argument( '-p', '--padPaired', dest='padPaired', action='store_true', default=False,
                        help='Pad paired reads so that they are the same length after trimming all trimming has occured. N\'s will be added to the 3\' end with \'#\' added to the quality score for each N that is added. This will not do anything for single-end reads. (default: no)' )
argParser.add_argument( '-m', '--minLen', dest='minLen', action='store', default=0, type=int,
                        help='Minimum size of trimmed read. If trimmed beyond minLen, then read is discarded. If read is paired then read is replaced with N\'s, unless both reads in pair are smaller than minLen in which case the pair is discarded. (default: no minimum length)' )
argParser.add_argument( '-c3', '--cut3', dest='num_cut_3', action='store', default='0', type=int, 
                        help='number of bases to remove from 3\' end of read. Truncating reads does not count toward trimming totals. (default: 0)' )
argParser.add_argument( '-c5', '--cut5', dest='num_cut_5', action='store', default='0', type=int, 
                        help='number of bases to remove from 5\' end of read. Truncating reads does not count toward trimming totals.  (default: 0)' )
argParser.add_argument( '-rN', '--removeNs', dest='removeN', action='store_true', default=False,
                        help='remove N\'s from both ends of the read. (default: no)' )
argParser.add_argument( '-c', dest='contaminants_fa', action='store', default=None, 
                        help='fasta-like file containing list of contaminants to trim' )
argParser.add_argument( '-f', dest='forward_fq', action='store', required=True,
                        help='fastq file with reads to be trimmed' )
argParser.add_argument( '-r', dest='reverse_fq', 
                        help='second read with paired-end reads. This file is not present for single-end reads.' )
argParser.add_argument( '-o', dest='output_prefix', action='store', required=True,
                        help='prefix for output file(s). A \'_1\' will be appended to the forward reads output file and if paired reads then a \'_2\' will be appended to the reverse reads file.' )

clArgs = argParser.parse_args()
if DEBUG: print clArgs

# flag if paired-end reads
PAIRED = False
if clArgs.reverse_fq:
    PAIRED = True

# track trimming stats
nPairsTrimmed = 0 # total number of trimmed read pairs (equal to numForTrimmed if single-end)
nForTrimmed = 0 # num forward reads trimmed
nRevTrimmed = 0 # num reverse reads trimmed
nPairsDiscarded = 0 # num read pairs discarded (equal to nForDiscarded if single-end)
nForDiscarded = 0 # num of forward reads replaced with N's
nRevDiscarded = 0 # num of reverse reads replaced with N's
nTotalReads = 0 # total number of reads

# print error and quit
def quitOnError(msg):
    print 'ERROR:', msg
    sys.exit(0)

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

    # trim N from beginning of read (5' end)
    if seq.startswith('N'):
        # trim sequence
        seq = seq.lstrip('N')

        # need to trim quals in same way we trimmed sequence
        quals = quals[len(quals) - len(seq):]

        # flag read as having been trimmed
        read[TRIMMED] = True

    # trim N from end of read (3' end)
    if seq.endswith('N'):
        # trim sequence
        seq = seq.rstrip('N')

        # need to trim quals in same way we trimmed sequence
        quals = quals[:len(seq)]

        # flag read as having been trimmed
        read[TRIMMED] = True

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    return read

# trim contaminant removing entire contaminant based on single k-mer
# mapping, from the inner end of the contaminant to the outer end of
# the read sequence.
def fullContaminantTrimming(read, contaminant):
    """
    Trim read based on contaminant. Contaminant is expected to be a
    tuple containing a set of trimming options and the contaminant
    sequence. For trimming first map a k-mer then trim the read from
    the expected start point of the contaminant based on the position
    of the k-mer within the contaminant to the end of the read. Note
    that this function is an agressive trimming scheme and may remove
    portions of the read that don't map to the contaminant.
    Contaminant example: [{'windows': 5, 'max': 0, 'end': 3, 'name':
    'oligo', 'size': 10}, 'AATTTAAATTTAATTTCCGGGGAATTANN',
    [('AATTTAAATT', 0), ('TAAATTTAAT', 4), ('TTTAATTTCC', 8),
    ('ATTTCCGGGG', 12), ('CCGGGGAATT', 16)]]
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
    # max number of bases to trim
    maxTrimmed = options['max']
    # which end to trim
    endTrimmed = options['end']

    # we attempt to perfectly align each k-mer against the
    # sequence. If any alignment succeeds, the index position of the
    # start or end of the original k-mer sequence is used for
    # trimming, depending on which end of the read sequence is
    # expected to contain the contaminant. For example if the mapped
    # k-mer begins at the 5th base in sequence A and it's expected to
    # map to the 3' end of sequence B then sequence A is assumed to
    # extend another 4 bases toward the 5' end of sequence B and hence
    # that position in sequence B is used for trimming. If none of the
    # k-mer's map then no trimming occurs.
    pos = -1
    for kmer, offset in kmers:
        # look for k-mer in read sequence. Index is the 5' position in
        # the read where the k-mer mapped
        index = seq.find(kmer)
        if index > -1: 
            # the read sequence may not contain the entire contaminant
            # so change the order of the k-mer search based on the end
            # we are searching. If we're expecting the k-mer to map to
            # the 3' end of the read then we search the k-mers that
            # are located near the 5' end of the contaminant before
            # searching for the k-mers located near the 3' end of the
            # contaminant.
            if endTrimmed == 3: pos = index - offset
            else: pos = index + offset

            # break out of loop since we found a match
            break

    # if k-mers not found, then no trimming
    if pos == -1: return read

    if endTrimmed == 3:
        seq = seq[:pos]
        quals = quals[:pos]
    else:
        seq = seq[pos:]
        quals = quals[pos:]

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    read[TRIMMED] = True
    return read

# trim contaminant removing the portion of the contaminant that maps
# from the k-mer region to the end of the sequence.
def mappedContaminantTrimming(read, contaminant):
    """
    Trim read based on contaminant. Contaminant is expected to be a
    tuple containing a set of trimming options and the contaminant
    sequence. For trimming first map a k-mer then extend the mapping
    toward the opposite end of the read. After the mapping is fully
    extended, trim the read from the contaminant mapping regaion to
    the end of the read. Note that this function is a more
    conservative trimming scheme and may miss contaminate regions that
    are smaller than the size of a k-mer. Contaminant example:
    [{'windows': 5, 'max': 0, 'end': 3, 'name': 'oligo', 'size': 10},
    'AATTTAAATTTAATTTCCGGGGAATTANN', [('AATTTAAATT', 0),
    ('TAAATTTAAT', 4), ('TTTAATTTCC', 8), ('ATTTCCGGGG', 12),
    ('CCGGGGAATT', 16)]]
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
    # max number of bases to trim
    maxTrimmed = options['max']
    # which end to trim
    endTrimmed = options['end']

    # we attempt to perfectly align each k-mer against the
    # sequence. Once we find a match, then we extend the mapping away
    # from the end of the read, only allowing for exact matches. We
    # then remove the mapped region of the read and any additional
    # bases on that end of the read. For example if the mapped k-mer
    # begins at the 5th base in sequence A and it's expected to map to
    # the 3' end of sequence B then we extend the k-mer mapped toward
    # the 5' end of the sequences.  The position returned is the last
    # mapped base on the 5' end of the mapping. If none of the k-mer's
    # map then no trimming occurs.
    pos = -1
    for kmer, offset in kmers:
        # when we set pos we break out of the inner loop. Here we
        # break out of the outer loop.
        if pos > -1: break

        # look for k-mer in read sequence. Index is the 5' position in
        # the read where the k-mer mapped
        index = seq.find(kmer)
        if index > -1: 
            # the read sequence may not contain the entire contaminant
            # so change the order of the k-mer search based on the end
            # we are searching. If we're expecting the k-mer to map to
            # the 3' end of the read then we search the k-mers that
            # are located near the 5' end of the contaminant before
            # searching for the k-mers located near the 3' end of the
            # contaminant.
            if endTrimmed == 3:
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
            else: 
                # found mapping so now extend mapping toward 3' end

                # need to adjust index since it's registering the 5'
                # end of the mapping and we are interested in the 3'
                # end of the mapping.
                index += options['size']

                pos = length
                while offset < cLength and index < length:
                    # offset is the position of the 3' end of the k-mer
                    # in the contaminant. So if equal to length of
                    # contaminant then we are already at the 3' end of
                    # the contaminant. index is our position in the
                    # read sequence. If that's equal to the length of
                    # the read sequence then we are at the 3' end of
                    # the read.
                    index += 1
                    offset += 1
                    if seq[index] != cSeq[offset]:
                        pos = index
                        break

    # if k-mers not found, then no trimming
    if pos == -1: return read

    if endTrimmed == 3:
        seq = seq[:pos]
        quals = quals[:pos]
    else:
        seq = seq[pos:]
        quals = quals[pos:]

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    read[TRIMMED] = True
    return read

# trim contaminant removing entire contaminant based on single k-mer
# mapping and whether or not the percent identity between the
# contaminant and the read is above the specified thresholds
def identityTrimming(read, contaminant):
    """
    Trim read based on contaminant. Contaminant is expected to be a
    tuple containing a set of trimming options and the contaminant
    sequence. For trimming first map a k-mer to align the contaminant
    to the read. Then compute the percent identity between the read
    and contaminant. If the percent identity is above the specified
    threshold then trim read. If the threshold is not reached then do
    not trim the read. With a high threshold this function will biast
    away from trimming. As the threshold is decreased (toward 0) then
    the likelyhood of trimming increases.  Contaminant example:
    [{'windows': 5, 'max': 0, 'end': 3, 'name': 'oligo', 'size': 10},
    'AATTTAAATTTAATTTCCGGGGAATTANN', [('AATTTAAATT', 0),
    ('TAAATTTAAT', 4), ('TTTAATTTCC', 8), ('ATTTCCGGGG', 12),
    ('CCGGGGAATT', 16)]]
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
    # max number of bases to trim
    maxTrimmed = options['max']
    # which end to trim
    endTrimmed = options['end']
    # percent identity threshold
    percIdentityThreshold = options['percentIdentity']
    # total identity threshold
    totIdentityThreshold = options['totalIdentity']

    # we attempt to perfectly align each k-mer against the
    # sequence. If any alignment succeeds, then we compute the percent
    # identity between the contaminant and read. If above our
    # threshold then trim, otherwise do not trim.
    percIdentity = 0
    totIdentity = 0
    for kmer, offset in kmers:
        # look for k-mer in read sequence. Index is the 5' position in
        # the read where the k-mer mapped
        index = seq.find(kmer)
        if index > -1: 
            # the read sequence may not contain the entire contaminant
            # so change the order of the k-mer search based on the end
            # we are searching. If we're expecting the k-mer to map to
            # the 3' end of the read then we search the k-mers that
            # are located near the 5' end of the contaminant before
            # searching for the k-mers located near the 3' end of the
            # contaminant.
            if endTrimmed == 3:
                # Pos is the position where the 5' end of contaminant
                # maps to read, based in k-mer mapping. This is the
                # position that will be trimmed if the percent
                # identity is high enough to allow for
                # trimming. Offset is the position of the 5' end of
                # the contaminant.
                pos = index - offset

                # cPos is the position in the contaminant
                cPos = 0
                # rPos is the read sequence position
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

            else: 
                # since we're now expecting the contaminant on the 5'
                # end, pos is the position where the 3' end of
                # contaminant maps to the read, based in k-mer
                # mapping. This is the position that will be trimmed
                # if the percent identity is high enough to allow for
                # trimming. Offset is the position of the 5' end of
                # the contaminant.
                pos = index + offset

                # cPos is the position in the contaminant
                cPos = cLength - (offset + index)
                # rPos is the read sequence position
                rPos = 0
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

    # if the k-mers were not found or the contaminant didn't map
    # sufficiently to the read then we don't trim.
    if (percIdentity < percIdentityThreshold) or (totIdentity < totIdentityThreshold): return read

    # if k-mers not found, then no trimming
    if pos == -1: return read

    if endTrimmed == 3:
        seq = seq[:pos]
        quals = quals[:pos]
    else:
        seq = seq[pos:]
        quals = quals[pos:]

    read[SEQUENCE] = seq
    read[QUALS] = quals
    read[LENGTH] = len(seq)
    read[TRIMMED] = True
    return read

#------------------------------------------------------------------------------------------
# LOAD CONTAMINANTS
#------------------------------------------------------------------------------------------

def computeKMers(seq, kmerSize, numWindows, end, method):
    """
    This will return a tuple of k-mers that includes the k-mer
    sequence. The tuple will also include either the 5' position of
    the k-mer within the original sequence or the 3' position of the
    k-mer within the original sequence, depending on which end the
    contaminant is expected to be located within the read
    sequence. Lastly, the order of the k-mers within the tuple will
    also depend on the read sequence end. We don't expect the read
    sequence to contain the entire contaminant so change the order of
    the k-mer search based on the end we are searching. If we're
    expecting the k-mer to map to the 3' end of the read then we
    search the k-mers that are located near the 5' end of the
    contaminant before searching for the k-mers located near the 3'
    end of the contaminant.
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
    windowStep = lastWindow / (numWindows - 1)

    kmerList = []
    # slice up the sequence based on the window step size and the k-mer size
    if end == 3:
        for index in range(0, lastWindow + 1, windowStep):
            # store k-mer and the k-mer position offset for a 3'
            # mapping. Offset is the 5' position of the k-mer.
            kmerList.append((seq[index:index+kmerSize], index))
    else:
        for index in range(lastWindow, 0, -1 * windowStep):
            if method == 0:
                # store k-mer and the k-mer position offset for a 5'
                # mapping using the Full trimming method. Offset is
                # the number of positions from the 5' end of the k-mer
                # to the 3' end of the contaminant.
                kmerList.append((seq[index:index+kmerSize], (seqLength-1) - (index-1)))
            elif method == 1:
                # store k-mer and the k-mer position offset for a 5'
                # mapping using the Mapped trimming method. Offset is
                # the 3' position of the k-mer.
                kmerList.append((seq[index:index+kmerSize], index+kmerSize))
            else:
                # store k-mer and the k-mer position offset for a 5'
                # mapping using the Percent Identity trimming
                # method. Offset is the 3' position of the k-mer.
                kmerList.append((seq[index:index+kmerSize], (seqLength-1) - (index-1)))
                #kmerList.append((seq[index:index+kmerSize], index+kmerSize))

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
            options = { 'end':DEFAULT_END, 
                        'max':DEFAULT_MAX, 
                        'method':DEFAULT_METHOD, 
                        'pecentIdentity':DEFAULT_PERCENT_IDENTITY, 
                        'size':DEFAULT_SIZE_KMER, 
                        'totalIdentity':DEFAULT_TOTAL_IDENTITY, 
                        'windows':DEFAULT_WINDOWS_NUMBER }
            for option in args:
                key, value = option.split(':')

                if key == 'end':
                    value = int(value)
                    if value != 3 and value != 5:
                        msg = 'Invalid "end" option (%d) for contaminant. Must be 3 or 5.' % value
                        quitOnError(msg)
                elif key == 'max':
                    value = int(value)
                elif key == 'method':
                    value = int(value)
                    if value != 0 and value != 1 and value != 2:
                        msg = 'Invalid "method" option (%d) for contaminant. Must be 0, 1, or 2.' % value
                        quitOnError(msg)
                elif key == 'percentIdentity':
                    value = float(value)
                    if value <= 0.0 or value > 1.0:
                        msg = 'Invalid "percentIdentity" option (%f) for contaminant. Percent identity must be greater than 0 and less than or equal to 1.' % value
                        quitOnError(msg)
                elif key == 'size':
                    value = int(value)
                elif key == 'totalIdentity':
                    value = int(value)
                elif key == 'windows':
                    value = int(value)
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

        else:
            seq = line.strip().upper()

            # compute set of k-mers
            kmerList = computeKMers(seq, options['size'], options['windows'], options['end'], options['method'])

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

forReadIn = ''
try: forReadIn = open(clArgs.forward_fq, 'r')
except: print 'ERROR: Unable to load reads file', clArgs.forward_fq

forReadOut = ''
try: forReadOut = open(clArgs.output_prefix + "_1.fq", 'w')
except: print 'ERROR: Unable to open output file', clArgs.output_prefix + "_1.fq"

if PAIRED:
    revReadIn = ''
    try: revReadIn = open(clArgs.reverse_fq, 'r')
    except: print 'ERROR: Unable to load reads file', clArgs.reverse_fq

    revReadOut = ''
    try: revReadOut = open(clArgs.output_prefix + "_2.fq", 'w')
    except: print 'ERROR: Unable to open output file', clArgs.output_prefix + "_2.fq"

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
            
    nTotalReads += 1

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
                forRead = fullContaminantTrimming(forRead, contaminant)
                if PAIRED: revRead = fullContaminantTrimming(revRead, contaminant)
            elif method == 1:
                forRead = mappedContaminantTrimming(forRead, contaminant)
                if PAIRED: revRead = mappedContaminantTrimming(revRead, contaminant)
            else:
                forRead = identityTrimming(forRead, contaminant)
                if PAIRED: revRead = identityTrimming(revRead, contaminant)

    #--------------------------------------------------------------------------------------
    # compute trimming stats
    trimFlag = False
    if forRead[TRIMMED]: 
        nForTrimmed += 1
        trimFlag = True
    if PAIRED: 
        if revRead[TRIMMED]: 
            nRevTrimmed += 1
            trimFlag = True
    # count paired reads as a single read in nPairsTrimmed
    if trimFlag: nPairsTrimmed += 1

    #--------------------------------------------------------------------------------------
    # discard reads that are too short
    if clArgs.minLen > 0:
        discardFor = False
        discardRev = False
        if forRead[LENGTH] < clArgs.minLen:
            nForDiscarded += 1
            discardFor = True
        if PAIRED:
            if revRead[LENGTH] < clArgs.minLen:
                nRevDiscarded += 1
                discardRev = True
        if discardFor or discardRev: nPairsDiscarded += 1

        # XXX could output them to an error file here
        if discardFor and discardRev:
            # both reads in read pair are discarded so don't write them to output file
            continue
        if discardFor and not PAIRED:
            # single-end read being discarded, so don't write to output file
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

print 'Reads processed:', nTotalReads
print 'Forward reads trimmed:', nForTrimmed
print 'Reverse reads trimmed:', nRevTrimmed
print 'Reads pairs trimmed:', nPairsTrimmed
print 'Forward reads discarded:', nForDiscarded
print 'Reverse reads discarded:', nRevDiscarded
print 'Read pairs discarded:', nPairsDiscarded

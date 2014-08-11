#!/usr/bin/env python

# Copyright (c) 2014, Stephen Fisher and Junhyong Kim, University of
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
by: S. Fisher, 2014

usage: removeDuplicates.py [-h] -f FIRST_FQ [-r SECOND_FQ] -o OUTPUT_PREFIX

This will return a file that contains the specified number of randomly
sampled lines from the original file. If 'lines grouped' is greater
than 1, then each time a line is selected, the specified number of
lines (grouping size) will also be include. For example a line
grouping of 4 means 4 lines will be included every time a line is
selected, as in the case of a FASTQ file.
"""

#------------------------------------------------------------------------------------------
# INITIALIZATIONS
#------------------------------------------------------------------------------------------

import sys, os, argparse

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

VERSION = '0.1'

# indecies for the read set
HEADER = 'header'
SEQUENCE = 'seq'
QUALS = 'quals'

argParser = argparse.ArgumentParser(version=VERSION, 
                                    description='Remove duplicate reads.',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    epilog='' +
                                    'Remove duplicate reads that exactly match.\n')
argParser.add_argument( '-f', dest='first_fq', action='store', required=True,
                        help='fastq file with reads to be processed' )
argParser.add_argument( '-r', dest='second_fq', 
                        help='fastq file with mate pairs for paired-end reads. This file is not present for single-end reads.' )
argParser.add_argument( '-o', dest='output_prefix', action='store', required=True,
                        help='prefix for output file(s). A \'_1\' will be appended to the first reads output file and if paired reads then a \'_2\' will be appended to the second reads file. Output files similarly named and with the suffix \'.dupl.txt\' will be created to store a copy of all duplicate reads.' )

clArgs = argParser.parse_args()
if DEBUG: print clArgs

# flag if paired-end reads
PAIRED = False
if clArgs.second_fq:
    PAIRED = True

# print error and quit
def quitOnError(msg):
    # XXX may not be compatible with Python v3
    sys.stderr.write('ERROR:' + msg + '\n')
    sys.exit(0)

# track stats
nInitialReadPairs = 0 # number of input read pairs
nUniqueReadPairs = 0 # number of output read pairs (ie after duplicates removed)
nDuplicatesFound = 0 # number of duplicate read pairs discarded
nReadsWithDuplicates = 0 # number of read pairs that contained at least 1 duplicate
maxDuplicateCount = 0 # max number of duplicates of a single read pair

uniqReads = {} # dictionary containing all uniq reads
duplReads = {} # dictionary containing a single copy of all duplicate reads
cntDupl = {} # dictionary containing a duplicates count for all reads

#------------------------------------------------------------------------------------------
# FUNCTIONS FOR HANDLING READING AND WRITING OF READS
#------------------------------------------------------------------------------------------

def nextRead(inFile):
    """
    Returns a set consisting of each line from the read. Returns empty set if no more reads.
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

        # Ignore the + in the fastq file
        inFile.readline()

        read[QUALS] = inFile.readline().rstrip()
        
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
# OPEN READ INPUT AND OUTPUT FILES
#------------------------------------------------------------------------------------------

# open input fastq file
firstReadIn = ''
try: 
    firstReadIn = open(clArgs.first_fq, 'r')
except: 
    msg = 'Unable to load reads file ' + clArgs.first_fq
    quitOnError(msg)

# open output fastq file
firstReadOut = ''
try: 
    firstReadOut = open(clArgs.output_prefix + "_1.fq", 'w')
except: 
    msg = 'Unable to open output file ' + clArgs.output_prefix + "_1.fq"
    quitOnError(msg)

# open file that will store one copy of all duplicate reads
firstReadDuplicateOut = ''
try: 
    firstReadDuplicateOut = open(clArgs.output_prefix + "_1.dupl.fq", 'w')
except: 
    msg = 'Unable to open output file ' + clArgs.output_prefix + "_1.dupl.fq"
    quitOnError(msg)

if PAIRED:
    secondReadIn = ''
    try: 
        secondReadIn = open(clArgs.second_fq, 'r')
    except: 
        msg = 'Unable to load reads file ' + clArgs.second_fq
        quitOnError(msg)

    secondReadOut = ''
    try: 
        secondReadOut = open(clArgs.output_prefix + "_2.fq", 'w')
    except: 
        msg = 'Unable to open output file ' + clArgs.output_prefix + "_2.fq"
        quitOnError(msg)

    secondReadDuplicateOut = ''
    try: 
        secondReadDuplicateOut = open(clArgs.output_prefix + "_2.dupl.fq", 'w')
    except: 
        msg = 'Unable to open output file ' + clArgs.output_prefix + "_2.dupl.fq"
        quitOnError(msg)

# generate tsv file with number of duplicates per read.
if DEBUG:
    readCountsOut = ''
    try:
        print clArgs.output_prefix + ".debug.tsv"
        readCountsOut = open(clArgs.output_prefix + ".debug.tsv", 'w')
    except: 
        msg = 'Unable to open output file ' + clArgs.output_prefix + '.debug.tsv'
        quitOnError(msg)

#------------------------------------------------------------------------------------------
# PROCESS READS. LOAD ONE READ FROM BOTH FQ FILES AT SAME
# TIME. PROCESS READ PAIR TOGETHER.
# ------------------------------------------------------------------------------------------

while 1:
    firstRead = nextRead(firstReadIn)
    if firstRead == {}: break
    readSeq = firstRead[SEQUENCE]
    if PAIRED: 
        secondRead = nextRead(secondReadIn)
        readSeq += secondRead[SEQUENCE]
    else:
        # create empty dict so don't have to keep testing PAIRED below
        secondRead = {}

    nInitialReadPairs += 1

    if readSeq in uniqReads:
        # found a duplicate
        if readSeq not in duplReads:
            # new duplicate, so add to duplRead and cntDupl
            duplReads[readSeq] = [firstRead, secondRead]
            cntDupl[readSeq] = 1
        else:
            # already found a duplicate so just increment cntDupl
            cntDupl[readSeq] += 1
    else:
        uniqReads[readSeq] = [firstRead, secondRead]

# we've completed the hash tables at this point so now we need to write the data
if PAIRED:
    for value in uniqReads.itervalues():
        writeRead(value[0], firstReadOut)
        writeRead(value[1], secondReadOut)

    for value in duplReads.itervalues():
        writeRead(value[0], firstReadDuplicateOut)
        writeRead(value[1], secondReadDuplicateOut)
else:
    for value in uniqReads.itervalues():
        writeRead(value[0], firstReadOut)

    for value in duplReads.itervalues():
        writeRead(value[0], firstReadDuplicateOut)

firstReadIn.close()
firstReadOut.close()
firstReadDuplicateOut.close()
if PAIRED:
    secondReadIn.close()
    secondReadOut.close()
    secondReadDuplicateOut.close()

if DEBUG:
    readCountsOut.write('numDupl\theader\tmate1\tmate2\n')
    for key in cntDupl:
        readCountsOut.write(str(cntDupl[key]) + '\t' + duplReads[key][0][HEADER] + '\t' + duplReads[key][0][SEQUENCE] + '\t' + duplReads[key][1][SEQUENCE] + '\n')
    readCountsOut.close()

# compute stats
nUniqueReadPairs = len(uniqReads)
nReadsWithDuplicates = len(duplReads)
nDuplicatesFound = sum(cntDupl.values())
if cntDupl != {}:
    maxDuplicateCount = max(cntDupl.values())

#------------------------------------------------------------------------------------------
# OUTPUT STATS
#------------------------------------------------------------------------------------------

print 'Reads processed:', nInitialReadPairs
print 'Unique reads:', nUniqueReadPairs
print 'Duplicate reads discarded:', nDuplicatesFound
print 'Number reads with at least 1 duplicate:', nReadsWithDuplicates
print 'Max duplicates of a single read:', maxDuplicateCount

# tab delimited output to facilitate adding stats to compilation file
fields = '\nnInitialReadPairs\tnUniqueReadPairs\tnDuplicatesFound\tnReadsWithDuplicates\tmaxDuplicateCount'
counts = '%d\t%d\t%d\t%d\t%d' % (nInitialReadPairs, nUniqueReadPairs, nDuplicatesFound, nReadsWithDuplicates, maxDuplicateCount)

print fields
print counts

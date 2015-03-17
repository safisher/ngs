#!/usr/bin/env python

# Copyright (c) 2015, Stephen Fisher and Junhyong Kim, University of
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
by: S. Fisher, 2015

usage: combineGeneCounts.py <AND | OR> <FILE 1> <FILE 2> <OUTPUT FILENAME>
"""

import sys, os

#------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

# expect 4 args
if len(sys.argv) < 5:
    print 'Usage: combineGeneCounts.py  <AND | OR> <FILE_1> <FILE_2> <OUTPUT FILE>'
    print '\tCombine two files containing gene counts according to the boolean operator specified.\n'
    print '\tThe boolean operator must be either "AND" or "OR".'
    print '\t\tAND: output count is zero if either file has a zero count.'
    print '\t\tOR: output count is equal to the sum of the counts from both files.\n'
    print '\tThe files are formatted as follows:'
    print '\t\tThe first row in each file is expected to be column headers (name count).'
    print '\t\tColumns are expected to be tab delimited.\n'
    print '\tAny genes in FILE_1 that are not in FILE_2 will be skipped. Any genes in FILE_2 '
    print '\tthat are not in FILE_1 will be treated as if they had a read count of "0" in FILE_1.\n'
    print '\tThe output file will be overwritten without warning.\n'
    print '\tExample usage:\n\t\tcombineGeneCounts.py AND exonCnts intronCnts combinedCnts'
    sys.exit()

BOOLEAN = sys.argv[1].lower()
FILE_1 = sys.argv[2]
FILE_2 = sys.argv[3]
OUT_FILE = sys.argv[4]

file1 = open(FILE_1, 'r')
file2 = open(FILE_2, 'r')
outFile = open(OUT_FILE, 'w')

# make sure boolean is either AND or OR
if BOOLEAN not in ('and', 'or'):
    print 'ERROR: Incorrect boolean operator: ' + BOOLEAN
    print 'Boolean operator must be either "AND" or "OR".'
    sys.exit()

# build dictionary from first file. We load the entire file into RAM,
# so this is limited by the available RAM.
geneCounts = {}
count = 0
firstLine = True
for line in file1:
    #if DEBUG: print line.rstrip()

    # skip header line
    if firstLine:
        firstLine = False
        continue

    # skip header line
    gene = line.split('\t')
    geneID = gene[0]
    geneCounts[geneID] = int(gene[1])  # read count
    count = count + 1

    if DEBUG: print count, geneID, gene[1]

# add header to output file
outFile.write('gene\tcount\n')

firstLine = True
for line in file2:
    # skip header line
    if firstLine:
        firstLine = False
        continue

    gene = line.split('\t')
    geneID = gene[0]
    readCnt_2 = int(gene[1])

    if DEBUG: print geneID, readCnt_2

    # make sure the gene exists in first file. If not then treat as if
    # it has a 0 value in the first file.
    if geneID not in geneCounts:
        print 'WARNING: Gene "' + geneID + '" is not in "' + FILE_1 + '". TREATING GENE AS IF WAS "0" IN THE FIRST FILE.'
        readCnt_1 = 0

    else:
        # get read count from first file
        readCnt_1 = geneCounts[geneID]

        # remove gene from list, so we can test if all genes from first
        # file are in the second file.
        geneCounts.pop(geneID)

    # combine the counts and output based on boolean operator
    finalCount = 0
    if BOOLEAN == 'and':
        # for AND both files must have a non-zero count for the output
        # to have a non-zero count.
        if readCnt_1 > 0 and readCnt_2 > 0:
            finalCount = readCnt_1 + readCnt_2
    else:
        finalCount = readCnt_1 + readCnt_2

    # output the combined count
    outFile.write(geneID + '\t' + str(finalCount) + '\n')

if geneCounts:
    # still have items in geneCounts so files not equal
    print 'ERROR: The following genes are not in the file "' + FILE_2 + '".'
    print 'THESE GENES ARE NOT INCLUDED IN THE OUTPUT FILE.'
    print geneCounts

file1.close()
file2.close()
outFile.close()

print 'Processed ' + str(count) + ' genes.'

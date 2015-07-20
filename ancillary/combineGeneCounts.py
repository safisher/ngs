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

usage: combineGeneCounts.py <AND | OR> <FILE 1> <FILE 2>
"""

import sys, os

#------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

# disable warnings
QUIET = True

# expect 3 args
if len(sys.argv) < 4:
    sys.stderr.write('Usage: combineGeneCounts.py  <AND | OR> <FILE_1> <FILE_2>\n')
    sys.stderr.write('Combine two files containing gene counts according to the boolean operator specified.\n\n')
    sys.stderr.write('The boolean operator must be either "AND" or "OR".\n')
    sys.stderr.write('\tAND: output count is zero if either file has a zero count.\n')
    sys.stderr.write('\t\tAny genes in FILE_1 that are not in FILE_2 will be skipped. Any genes in FILE_2 \n')
    sys.stderr.write('\t\tthat are not in FILE_1 will be treated as if they had a read count of "0" in FILE_1.\n\n')
    sys.stderr.write('\tOR: output count is equal to the sum of the counts from both files.\n\n')
    sys.stderr.write('\t\tIf a gene is missing from one of the files, it is assumed to have a "0" value in that file\n\n')
    sys.stderr.write('The files are formatted as follows:\n')
    sys.stderr.write('\tThe first row in each file is expected to be column headers (name count).\n')
    sys.stderr.write('\tColumns are expected to be tab delimited.\n\n')
    sys.stderr.write('The output will be written to the screen (i.e. STDOUT).\n\n')
    sys.stderr.write('Example usage:\n\t\tcombineGeneCounts.py OR exonCnts intronCnts\n')
    sys.exit()

BOOLEAN = sys.argv[1].lower()
FILE_1 = sys.argv[2]
FILE_2 = sys.argv[3]

if FILE_1 == '-':
    file1 = sys.stdin
else:
    file1 = open(FILE_1, 'r')

if FILE_2 == '-':
    file2 = sys.stdin
else:
    file2 = open(FILE_2, 'r')

# make sure boolean is either AND or OR
if BOOLEAN not in ('and', 'or'):
    sys.stderr.write('ERROR: Incorrect boolean operator: ' + BOOLEAN + '\n')
    sys.stderr.write('Boolean operator must be either "AND" or "OR".\n')
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

# send header to standard output
sys.stdout.write('gene\tcount\n')

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
        if not QUIET:
            sys.stderr.write('WARNING: Gene "' + geneID + '" is not in "' + FILE_1 + '". TREATING GENE AS IF WAS "0" IN THE FIRST FILE.\n')
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
    sys.stdout.write(geneID + '\t' + str(finalCount) + '\n')

# check to see if any genes were in file_1 and not file_2. If AND then
# this is an error. If OR then we assume they have 0 count in file_2
# and we'll output them.
if geneCounts:
    if BOOLEAN == 'and':
        # still have items in geneCounts so files not equal
        sys.stderr.write('ERROR: The following genes are not in the file "' + FILE_2 + '".\n')
        sys.stderr.write('THESE GENES ARE NOT INCLUDED IN THE OUTPUT FILE.\n')
        for gene in geneCounts:
            sys.stderr.write(gene)
        sys.stderr.write('\n')
    else:
        for gene in geneCounts:
            if not QUIET:
                sys.stderr.write('WARNING: Gene "' + gene + '" is not in "' + FILE_2 + '". TREATING GENE AS IF WAS "0" IN THE SECOND FILE.\n')
            sys.stdout.write(gene + '\t' + str(geneCounts[gene]) + '\n')

file1.close()
file2.close()

if DEBUG: print 'Processed ' + str(count) + ' genes.'

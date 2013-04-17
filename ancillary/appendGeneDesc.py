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

usage: appendGeneDesc.py <COUNTS> <ANNOTATION LIST> <OUTPUT PREFIX>
"""

import sys, os

#------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------

DEBUG = 0
if DEBUG: print 'DEBUG MODE: ON'

# expect 3 args
if len(sys.argv) < 4:
    print 'Usage: appendGeneDesc.py <COUNTS FILE> <ANNOTATION LIST> <OUTPUT FILE>'
    print '\tThis will append gene name and description to the *.cnts file that'
    print '\tis output from HTSEQ or runDESeq.py. The annotation list is expected'
    print '\tto be tab delimited with the following field order:'
    print '\t\tGene ID, Name, Description'
    sys.exit()

CNTS_FILE = sys.argv[1]
AN_FILE = sys.argv[2]
OUT_FILE = sys.argv[3]

cntsFile = open(CNTS_FILE, 'r')
anFile = open(AN_FILE, 'r')
outFile = open(OUT_FILE, 'w')

# build dictionary of annotations. We load the entire ANNOTATION LIST
# into RAM, so this is limited by the available RAM.
annotations = {}
count = 0
for line in anFile:
    #if DEBUG: print line.rstrip()

    annot = line.split('\t')
    geneID = annot[0]
    annotations[geneID] = [annot[1], annot[2]]  # gene name, gene desc
    count = count + 1

    #if DEBUG: print count, geneID, desc

print 'Processed ' + str(count) + ' annotations.'

count = 0
for line in cntsFile:
    if count == 0:
        outFile.write(line.strip() + '\tName\tDescription\n')
        count = 1
        continue

    geneID = line.split('\t')[0]

    # remove quotes around name
    #geneID = geneID[1:-1]
    geneID = geneID.replace('"', '').strip()

    # remove new line
    geneDesc = annotations[geneID][1].strip()

    # print to file
    outFile.write(line.strip() + '\t' + annotations[geneID][0] + '\t' + geneDesc + '\n')
    count = count + 1

print 'Processed ' + str(count) + ' genes.'

cntsFile.close()
anFile.close()
outFile.close()


#!/usr/bin/python

# Copyright (c) 2012, Stephen Fisher and Junhyong Kim, University of
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

usage: flyBase2HTSeq.py <exons GFF> <FBgn_FBtr> <Output>

Convert FlyBase GFF file to gene model appropriate for HTSeq
"""

import sys, os

#------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

# expect 3 args
if len(sys.argv) < 4:
    print 'Usage: flyBase2HTSeq.py <exons GFF> <FBgn_FBtr> <Output>'
    print '\texons GFF - list of exons from dmel-all-no-analysis-r5.49.gff'
    print '\tFBgn_FBtr - table gene and transcript IDs (fbgn_fbtr_fbpp_fb_2013_01.tsv)'
    print '\toutput - output file'
    sys.exit()

EXONS_FILE = sys.argv[1]
GN_TR_FILE = sys.argv[2]
OUT_FILE = sys.argv[3]

exonFile = open(EXONS_FILE, 'r')
gn2trFile = open(GN_TR_FILE, 'r')
outFile = open(OUT_FILE, 'w')

# build dictionary from gn2trFile. We load the entire list into RAM, so
# this is limited by the available RAM. File is expected to be two tab-delimited colums: geneID transcriptID
transcripts = {}
if DEBUG: tcounts = {}  # count every time we find an exon in our GFF file
count = 0
for line in gn2trFile:
    gene = line.split('\t')
    transcriptID = gene[1].strip() # need to remove \n from name
    transcripts[transcriptID] = gene[0] # gene ID
    if DEBUG: tcounts[transcriptID] = 0 # reset all counts
    count = count + 1

print 'Processed ' + str(count) + ' transcripts.'

count = 0
for line in exonFile:
    # GFF file is tab delimited
    exon = line.split('\t')
    attribs = exon[8].split(';') # semicolon separated list of attributes (eg ID=FBgn0031208:2;Parent=FBtr0300690)
    
    # get transcript ID(s) from attributes. Each exon might be in multiple transcripts
    # ex: ID=FBgn0003963:13;Name=ush:13;Parent=FBtr0078063,FBtr0329895;parent_type=mRNA
    for attrib in attribs:
        # the Parent key is the transcript ID. Skip other attributes
        if 'Parent=' not in attrib: continue

        # get list of transcript IDs
        # ex: Parent=FBtr0078063,FBtr0329895
        ids = (attrib.split('='))[1].split(',')
    
        for name in ids:
            # output exon to GFF file
            if name not in transcripts:
                print 'ERROR: Transcript ID (' + name + ') has no associated gene ID:\n' + line
                sys.exit()
            else:
                if DEBUG: tcounts[name] = tcounts[name] + 1
                # include transcript ID
                outFile.write(exon[0] + '\t' + exon[1] + '\t' + exon[2] + '\t' + exon[3] + '\t' + exon[4] + '\t' + exon[5] + '\t' + exon[6] + '\t' + exon[7] + '\tgene_id=' + transcripts[name] + ';trans_id=' + name + '\n')

    count = count + 1

print 'Processed ' + str(count) + ' exons.'

if DEBUG:
    for name in tcounts.keys():
        print transcripts[name] + '\t' + name + '\t' + str(tcounts[name])

exonFile.close()
gn2trFile.close()
outFile.close()


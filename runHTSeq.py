#!/usr/bin/env python

# Copyright (c) 2013, Hannah Dueck, Stephen Fisher and Junhyong Kim,
# University of Pennsylvania.  All Rights Reserved.
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
by: H. Dueck, S. Fisher, 2013

usage: python runHTSeq.py <SOURCE BAM> <OUTPUT FILE PREFIX> <GENE MODEL GTF>

Requires HTSeq version 0.6 or later as runHTSeq.py uses the --order flag and a position sorted SAM file.

HTSeq parameters: --mode=intersection-nonempty --stranded=no --type=exon --idattr=gene_id

Example:
  runHTSeq.py RUM_Unique.bam Sample_xeno3 /lab/repo/resources/htseq/zebrafish.gz

The goal here is to quantify reads in preprocessed BAM files 
(unique reads only, with paired end names identical, sorted by 
readname). 
The individual steps are to: 
  1) Quantify (X)
  2) synthesize quantified output.
"""

import sys, os, subprocess, datetime, pysam

#------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

# expect 3 args
if len(sys.argv) < 4:
    print 'Usage: python runHTSeq.py <SOURCE BAM> <OUTPUT FILE PREFIX> <GENE MODEL GTF>\n'
    sys.exit()

SOURCE_BAM = sys.argv[1]
SAMPLE = sys.argv[2] + '.htseq'
GENE_MODEL = sys.argv[3]

def count_features(sortedSam):
    # saves a tab separated file with feature names as row names, and overlap counts as values.
    # Note that this currently uses the intersection-nonempty method, counts only hits to exons and assigns counts to gene_ids.
    if DEBUG: update('Counting.')

    try:
	##MK: added the option --order=pos which runs HTseq on position sorted file
        counts = run_cmd('python -m HTSeq.scripts.count --order=pos --mode=intersection-nonempty --stranded=no --type=exon --idattr=gene_id ' + sortedSam + ' ' + GENE_MODEL)
    except:
        print 'Error during counting process: '
        print sys.exc_info()
        # don't delete the SAM file as that might help in debugging.
        update('Failed_counting')
        return 
	
    countsFile = SAMPLE + '.out'
    outFile = open(countsFile, 'w')
    outFile.writelines(counts)
    outFile.close()
	
    if DEBUG: update('Completed counting: ' + countsFile)

def update(message):
    print message
    print datetime.datetime.now()

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    return p.communicate()[0]

try:
    # convert source BAM to SAM
    sortedSam = SAMPLE + '.sorted.sam'
    run_cmd('samtools view -h -o ' + sortedSam + ' ' + SOURCE_BAM)
    count_features(sortedSam)

    # remove temporary file
    if not DEBUG: os.remove(sortedSam)
    
    if DEBUG: update('Finished.')
except:
    update('Error during processing: ')
    print  sys.exc_info()


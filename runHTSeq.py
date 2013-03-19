#!/usr/local/bin/python

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

usage: runHTSeq.py <SOURCE BAM> <OUTPUT FILE PREFIX> <GENE MODEL GTF>

HTSeq requires a SAM file that is sorted by read names and had reads fixed with samtools 'fixmate'. This script handles the conversion from an unsorted BAM file.

Example:
  runHTSeq.py RUM_Unique.bam Sample_xeno3 /lab/repo/resources/htseq/zebrafish/zebrafish.Zv9.69.gtf.gz

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

DEBUG = 0

if DEBUG: print 'DEBUG MODE: ON'

# expect 3 args
if len(sys.argv) < 4:
    print 'Usage: runHTSeq.py <SOURCE BAM> <OUTPUT FILE PREFIX> <GENE MODEL GTF>\n'
    sys.exit()

SOURCE_BAM = sys.argv[1]
SAMPLE = sys.argv[2] + '.htseq'
GENE_MODEL = sys.argv[3]

def count_features(sortedSam):
    # saves a tab separated file with feature names as row names, and overlap counts as values.
    # Note that this currently uses the intersection-nonempty method, counts only hits to exons and assigns counts to gene_ids.
    if DEBUG: update('Counting.')

    try:
        counts = run_cmd('python -m HTSeq.scripts.count --mode=intersection-nonempty --stranded=no --type=exon --idattr=gene_id ' + sortedSam + ' ' + GENE_MODEL)
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

def rename_reads():
    # Remove the suffix that RUM adds to read names (_a, _b, a or b).
    # Could add check to this method to skip unneeded loops.
    if DEBUG: update('Renaming pairs.')

    try:
        renamedBam = SAMPLE + '.renamed.bam'
        samfile = pysam.Samfile(SOURCE_BAM, "rb")
        renamed_reads = pysam.Samfile(renamedBam, "wb", template=samfile)
		
        for read in samfile.fetch():
            if read.qname[len(read.qname) - 2: len(read.qname)] == '_a':
                read.qname = read.qname[0:len(read.qname) - 2]
            elif read.qname[len(read.qname) - 2: len(read.qname)] == '_b':
                read.qname = read.qname[0:len(read.qname) - 2]
            elif read.qname[len(read.qname) - 1: len(read.qname)] == 'a':
                read.qname = read.qname[0:len(read.qname) - 1]
            elif read.qname[len(read.qname) - 1: len(read.qname)] == 'b':
                read.qname = read.qname[0:len(read.qname) - 1]
        renamed_reads.write(read)	
		
        renamed_reads.close()
        samfile.close()
		
    except:
        print 'Error during renaming of reads: '
        print sys.exc_info()

    if DEBUG: update('Completed renaming: ' + renamedBam)

    return renamedBam


def sort_by_readname(bamToSort):
    # Given a path to a bam file, saves a bam file sorted by read name as outpath_name_sorted.bam.
    if DEBUG: update('Sorting by read name.')

    # sort by read name. note that samtools will add ".bam" to output file name
    sortedBam = SAMPLE + '.sorted'
    run_cmd('samtools sort -n ' + bamToSort + ' ' + sortedBam)

    # samtools adds '.bam' to the file name
    sortedBam = SAMPLE + '.sorted.bam'
    
    if DEBUG: update('Completed sorting: ' + sortedBam)
    
    return sortedBam

def update(message):
    print message
    print datetime.datetime.now()

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    return p.communicate()[0]

try:
    # remove RUM suffix
    renamedBam = rename_reads()
    
    # sort BAM by read name
    sortedBam = sort_by_readname(renamedBam)

    # run fixmate on sorted BAM
    fixedBam = SAMPLE + '.fixed.bam'
    run_cmd('samtools fixmate ' + sortedBam + ' ' + fixedBam)

    # convert sorted and fixed BAM to SAM
    sortedSam = SAMPLE + '.sorted.sam'
    run_cmd('samtools view -h -o ' + sortedSam + ' ' + fixedBam)

    count_features(sortedSam)

    # remove temporary files
    os.remove(renamedBam)
    os.remove(sortedBam)
    os.remove(fixedBam)
    os.remove(sortedSam)
    
    if DEBUG: update('Finished.')
except:
    update('Error during processing: ')
    print  sys.exc_info()


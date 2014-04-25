#!/usr/bin/env python

# Copyright (c) 2012,2013, Stephen Fisher and Junhyong Kim, University of
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
by: S. Fisher, 2012

usage: python randomSample.py <num lines> <lines grouped> <input> <output>

This will return a file that contains the specified number of randomly
sampled lines from the original file. If 'lines grouped' is greater
than 1, then each time a line is selected, the specified number of
lines (grouping size) will also be include. For example a line
grouping of 4 means 4 lines will be included every time a line is
selected, as in the case of a FASTQ file.
"""

import sys, os, subprocess, random

#------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------

DEBUG = 0

if DEBUG: print 'DEBUG MODE: ON'

# expect 2 args
if len(sys.argv) < 3:
    print 'Usage: python randomSample.py <num lines> <lines grouped> <input> <output>'
    sys.exit()


NUM_LINES = int(sys.argv[1])
LINES_GROUPED = int(sys.argv[2])
IN_FILE = sys.argv[3]
OUT_FILE = sys.argv[4]

inFile = open(IN_FILE, 'r')
outFile = open(OUT_FILE, 'w')

# get number of lines
def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0: raise IOError(err)
    return int(result.strip().split()[0])

def error():
    print "Ran out of lines in the file\n"
    inFile.close()
    outFile.close()
    quit()

totalLines = file_len(IN_FILE) / LINES_GROUPED
lines = random.sample(xrange(totalLines), NUM_LINES)
lines.sort()

if DEBUG: print "totalLines: ", totalLines, "NUM_LINES: ", NUM_LINES
if DEBUG: print "lines: ", lines

i = 0
j = 0
numSampled = 0
for lnum in lines:
    while i < lnum:
        # assume lines are grouped, as in FASTQ files.
        g = 0
        while g < LINES_GROUPED:
            if DEBUG: print "searching j: ", j, "i: ", i, "g: ", g
            line = inFile.readline()
            if not line: error()
            g += 1
            j += 1
        i += 1

    if DEBUG: print "found one", i
    g = 0
    while g < LINES_GROUPED:
        if DEBUG: print "found j: ", j, "i: ", i, "g: ", g
        line = inFile.readline()
        if not line: error()
        outFile.write(line)
        g += 1
        j += 1
    i += 1
        
    numSampled += 1 

inFile.close()
outFile.close()

print "Total reads:", totalLines
print "Reads Sampled:", numSampled

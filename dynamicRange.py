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

usage: dynamicRange.py [-h] [-min min] [-max max] -c counts

This compute the dynamic range of a set of gene counts.
"""

#------------------------------------------------------------------------------------------
# INITIALIZATIONS
#------------------------------------------------------------------------------------------

import sys, os, argparse, math

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

VERSION = '0.1'

# indecies for the read set
HEADER = 'header'
SEQUENCE = 'seq'
QUALS = 'quals'

argParser = argparse.ArgumentParser(version=VERSION, 
                                    description='Compute dynamic range.',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    epilog='' +
                                    'This compute the dynamic range of a set of gene counts: log10((min/max)/(sum of all counts)). If max range is greater than the total number of genes, then the last gene count will be used (ie max=1500 but only 1200 genes so count for genea at position 1200 is used).\n')
argParser.add_argument( '-min', dest='rangeMin', default="150", action='store', 
                        help='lower bound of range' )
argParser.add_argument( '-max', dest='rangeMax', default="1500", action='store', 
                        help='upper bound of range' )
argParser.add_argument( '-c', dest='counts', required=True,
                        help='file containing gene counts. Expect two tab delimited columns (genes counts).' )

clArgs = argParser.parse_args()
if DEBUG: print clArgs

# subtract 1 because list is 0 based
rangeMin = int(clArgs.rangeMin) - 1
rangeMax = int(clArgs.rangeMax) - 1

#------------------------------------------------------------------------------------------
# OPEN INPUT FILE
#------------------------------------------------------------------------------------------

# open input counts file
countsIn = ''
try: 
    countsIn = open(clArgs.counts, 'r')
except: 
    msg = 'Unable to load counts file ' + clArgs.counts
    sys.stderr.write('ERROR:' + msg + '\n')
    sys.exit(0)

#------------------------------------------------------------------------------------------
# PROCESS READS. LOAD ONE READ FROM BOTH FQ FILES AT SAME
# TIME. PROCESS READ PAIR TOGETHER.
# ------------------------------------------------------------------------------------------

# expect first line to be "gene count" header
line = countsIn.readline().rstrip()
if len(line) == 0:
    sys.stderr.write('ERROR: Empty File\n')
    sys.exit(0)

counts = []
while 1:
    line = countsIn.readline().rstrip()
    if len(line) == 0:
        break

    cnt = int(line.split('\t')[1])

    # only include non-zero counts
    if cnt > 0:
        counts.append(cnt)

#------------------------------------------------------------------------------------------
# OUTPUT STATS
#------------------------------------------------------------------------------------------

counts.sort(reverse=True)
totalCounts = len(counts)

if totalCounts < rangeMax:
    dynamicRange = -1
else:
    dynamicRange = math.log(float(counts[rangeMin])/float(counts[rangeMax]))

if DEBUG:        
    if totalCounts < rangeMax:
        sys.stderr.write('ERROR: Only ' + str(totalCounts) + ' non-zero counts in file.\n')
        sys.exit(0)

    sumCounts = sum(counts)
    print 'Number of genes:', totalCounts
    print 'Sum of gene counts:', sumCounts
    print 'Rank 1 count:', counts[0]
    print 'Rank ' + str(rangeMin+1) + ' count:', counts[rangeMin]
    print 'Rank ' + str(rangeMax+1) + ' count:', counts[rangeMax]
    print 'Rank ' + str(totalCounts) + ' count:', counts[totalCounts-1]
    print 'Log Ratio:', math.log(float(counts[rangeMin])/float(counts[rangeMax]))

# tab delimited output to facilitate adding stats to compilation file
#fields = '\ndynamicRange'
#counts = '%f' % (dynamicRange)

#print fields
#print counts

print dynamicRange

#!/usr/bin/env python

# Copyright (c) 2011-2014, Stephen Fisher and Junhyong Kim, University of
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
by: S. Fisher, 2011,2014

Usage: python parseBlast.py targetSpecies readsFastaFile blastFile
"""

#------------------------------------------------------------------------------------------
# Initializations
#------------------------------------------------------------------------------------------

import sys, os

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

VERSION = '1.1'

# expect 3 args
if len(sys.argv) < 4:
    print 'Usage: python parseBlast.py targetSpecies readsFastaFile blastFile'
    print '\ttargetSpecies - reads mapped to this species will not be counted for other species. Target species must be one of the following species: bact, fish, fly, human, mouse, rat, yeast, ERCC'
    print '\tblast.csv - list of each hit per read and count of total number of hits per read'
    print '\treads.counted.txt - all alignments for reads that mapped to at least one of the target species'
    print '\treads.notCounted.txt - all alignments for reads that did not map to any target species'
    print '\tfailed.fa - fasta file with reads that did not align'
    print '\tfailed.tsv - blast results for read that did not map to any of the counted species'
    print '\tspeciesCounts.txt - species counts'
    print '\ttargetSpecies.tsv - blast results that mapped to target species'
    print '\ttargetSpecies.fa - fasta file with all target species reads'
    sys.exit()


TARGET = sys.argv[1]
RAW_FILE = sys.argv[2]
BLAST_FILE = sys.argv[3]
BLAST_PATH = os.path.dirname(RAW_FILE)
if len(BLAST_PATH) > 0:
    BLAST_PATH = BLAST_PATH + '/'

rawFile = open(RAW_FILE, 'r')
blastFile = open(BLAST_FILE, 'r')
csvFile = open(BLAST_PATH+'blast.csv', 'w')
hitFile = open(BLAST_PATH+'reads.counted.txt', 'w')
missFile = open(BLAST_PATH+'reads.notCounted.txt', 'w')
failedFile = open(BLAST_PATH+'failed.fa', 'w')
uncountedFile = open(BLAST_PATH+'uncounted.tsv', 'w')
speciesFile = open(BLAST_PATH+'speciesCounts.txt', 'w')

# convert from repo species name to names used here.
targetSpecies = TARGET.lower()
if 'mm9' in TARGET or 'mm10' in TARGET: targetSpecies = 'mouse'
if 'rn5' in TARGET: targetSpecies = 'rat'
if 'drosophila' in TARGET or 'dmel5' in TARGET: targetSpecies = 'fly'
if 'hg19' in TARGET: targetSpecies = 'human'
if 'hg38' in TARGET: targetSpecies = 'human'
if 'saccer3' in TARGET: targetSpecies = 'yeast'
if 'zebrafish' in TARGET or 'zv9' in TARGET: targetSpecies = 'fish'
if 'ercc' in TARGET: targetSpecies = 'ercc'

targetFile = open(BLAST_PATH+targetSpecies+'.tsv', 'w')
targetFileFa = open(BLAST_PATH+targetSpecies+'.fa', 'w')

#------------------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------------------

def error(line):
    # this should never happen
    print 'Got to this line in error: ', line, '\n'
    rawFile.close()
    blastFile.close()
    csvFile.close()
    hitFile.close()
    missFile.close()
    failedFile.close()
    uncountedFile.close()
    speciesFile.close()
    targetFile.close()
    targetFileFa.close()
    sys.exit()

def getLine():
    line = blastFile.readline()
    if not line: error(line)
    return line.rstrip()

def getRaw():
    # do this twice because first read is header info and we want sequence
    line0 = rawFile.readline()
    if not line0: error(line0)
    line1 = rawFile.readline()
    if not line1: error(line1)
    return line1.rstrip(), line0.rstrip()

def writeAlignments(outputFile):
    # output alignment info to outputFile
    outputFile.write('------------------------------------------------------------------------------------\n' )
    outputFile.write(query)
    for hit in range(0, hits):
        # step through until we find an alignment
        while True:
            line = getLine()

            # stop when we get to the alignment
            if '>' in line: break

            # if we don't find any alignments, then there's a problem with the blast file
            if 'Effective search space used' in line: error(line)

        # get alignement descriptor. It might be multiple lines and
        # will always begin with a '>' and end with a line with the length
        desc = ''
        while not 'Length=' in line:
            desc += line + ' '
            line = getLine()
        if len(desc) > 80: desc = desc[:77] + '...'
        outputFile.write('\n   ' + desc + '\n')

        # skip blank line
        getLine()

        # print line with E value
        line = getLine()
        outputFile.write('  ' + line + '\n')

        # print line with Identities and Gaps
        line = getLine()
        line = line.replace('(', '')
        line = line.replace('%)', ' ')
        line = line.replace('/', ' ')
        vals = line.split()

        # skip next 2 lines
        for cnt in range(0, 2): getLine()

        # an alignment has 3 lines
        for cnt in range(0, 3): 
            line = getLine()
            outputFile.write('      ' + line + '\n')
        
        # skip blank line
        line = getLine()

        # test if the alignment continues
        line = getLine()
        if 'Query' in line:
            outputFile.write('\n')
            for cnt in range(0, 3): 
                outputFile.write('      ' + line + '\n')
                line = getLine()

    outputFile.write('\n')

#------------------------------------------------------------------------------------------
# Do counts
#------------------------------------------------------------------------------------------

numReads = 0
numFails = 0
numHits = 0
numNotCounted = 0

# we count each time a read is aligned to one of the following. Since
# we only want to count the species once per read, we use a flag to
# denote whether or not the read has already been counted.
counts = {}
counts['bact'] = 0
counts['fish'] = 0
counts['fly'] = 0
counts['human'] = 0
counts['mouse'] = 0
counts['rat'] = 0
counts['yeast'] = 0
counts['ercc'] = 0

while True:
    line = blastFile.readline()
    if not line: break

    if not 'Query=' in line:
        # keep going until fine a 'Query='
        continue
    
    ###############################################
    # starting a new read
    query = line
    if DEBUG: print '\n--------------------------------------------------------------------------\n', line,

    rawSeq, seqHeader = getRaw()

    # count the read
    numReads += 1

    # skip next 5 lines
    for cnt in range(0, 5): line = getLine()

    # test if no alignments 
    if 'No hits' in line:
        if DEBUG: print '  Hits: 0'

        csvFile.write('0, ' + rawSeq + ', ' + seqHeader + '\n')

        # generate fasta file with sequences that didn't align
        failedFile.write(seqHeader + '\n' + rawSeq + '\n')
        
        # didn't align, so count the no hit and continue to next read
        numFails += 1

        # this will skip to next 'Query=' statement
        continue
    
    ###############################################
    # process hits

    # count number of hits for the read
    hits = 0

    # flag if we didn't count the read mapping
    hitNotCounted = False

    # when species is found then save the line so we can output to
    # targetFile (in the case of the target species). This also
    # servers to flag when a species has already been counted, so we
    # only count a species once per read.
    found = {}

    # flag if none of the hits were our species of interest
    notFound = ''

    while True:
        line = getLine()

        # if empty line then ran out of alignments
        if line == '': break

        hits += 1

        lLine = line.lower() # test lowercase version of line

        # keep first mapping as that will be the best mapping
        if ('bacter' in lLine) or ('coli' in lLine): 
            if 'bact' not in found: found['bact'] = line
        elif ('zebrafish' in lLine) or ('danio' in lLine) or ('rerio' in lLine):
            if 'fish' not in found: found['fish'] = line
        elif ('drosophila' in lLine) or ('melanogaster' in lLine):
            if 'fly' not in found: found['fly'] = line
        elif ('human' in lLine) or ('sapiens' in lLine): 
            if 'human' not in found: found['human'] = line
        elif ('mouse' in lLine) or ('musculus' in lLine): 
            if 'mouse' not in found: found['mouse'] = line
        elif ('rattus' in lLine) or ('norvegicus' in lLine): 
            if 'rat' not in found: found['rat'] = line
        elif ('yeast' in lLine) or ('cerevisiae' in lLine) or ('carlsbergensis' in lLine):
            if 'yeast' not in found: found['yeast'] = line
        elif ('ercc' in lLine):
            if 'ercc' not in found: found['ercc'] = line
        else:
            # count how many reads only mapped to species we are not tracking
            if len(notFound) == 0: notFound = line

        # print alignment
        if DEBUG: print ' ', line

    # finished processing read's hits so now figure out if target was mapped
    if targetSpecies in found:
        counts[targetSpecies] += 1
        
        targetFile.write(seqHeader + '\t' + rawSeq + '\t' + found[targetSpecies] + '\n')
        targetFileFa.write(seqHeader + '\n' + rawSeq + '\n')
    else:
        if found:
            for species in found.keys():
                counts[species] += 1
        else:
            numNotCounted += 1
            hitNotCounted = True
            uncountedFile.write(line + '\t' + rawSeq + '\t' + seqHeader + '\n')
            
    if DEBUG: print '  Hits:', hits
    csvFile.write(str(hits) + ', ' + rawSeq + ', ' + seqHeader + '\n')

    # count number of hits
    if hits > 0: numHits += 1

    ###############################################
    # process multiple alignments

    # XXX If this is set to '<2' then only multiple alignments will be included in hits file
    # if no alignments, then start on next query
    if hits == 0: continue
    
    if hitNotCounted: 
        writeAlignments(missFile)
    else: 
        writeAlignments(hitFile)


#------------------------------------------------------------------------------------------
# Output counts to species file
#------------------------------------------------------------------------------------------

speciesFile.write('parseBlast.py version: ' + VERSION + '\n')
speciesFile.write('Target Species: ' + targetSpecies + '\n')
speciesFile.write('Num reads: ' + str(numReads) + '\n')
speciesFile.write('Num reads that did not align: ' + str(numFails) + '\n')
speciesFile.write('Number hits: %d ( %.1f%% )\n' % (numHits, (100*float(numHits)/float(numReads))))
speciesFile.write('Hits not accounted for: ' + str(numNotCounted) + '\n')

# need to include ERCC in target counts when computing hits not target species.
targetCounts = counts[targetSpecies]
if 'ercc' not in 'targetSpecies': targetCounts += counts['ercc']
hitsNotTargetOrERCC = 100.0 * float(numHits - targetCounts)/float(numHits)
speciesFile.write('Hits not target species or ERCC: %.1f%%\n\n' % (hitsNotTargetOrERCC))

# print out species counts
keyList = sorted(counts.keys())
for species in keyList:
    speciesFile.write(species + ' hits\t\t' + str(counts[species]) + '\n')

speciesFile.write('\nTotal Hits\tHits Not Counted\tHits Not Target or ERCC\tBacteria\tFish\tFly\tHuman\tMouse\tRat\tYeast\tERCC\n')
speciesFile.write(str(numHits) + '\t' + str(numNotCounted))
speciesFile.write('\t%.1f%%' % (hitsNotTargetOrERCC)
speciesFile.write('\t' + str(counts['bact']))
speciesFile.write('\t' + str(counts['fish']))
speciesFile.write('\t' + str(counts['fly']))
speciesFile.write('\t' + str(counts['human']))
speciesFile.write('\t' + str(counts['mouse']))
speciesFile.write('\t' + str(counts['rat']))
speciesFile.write('\t' + str(counts['yeast']))
speciesFile.write('\t' + str(counts['ercc']))
speciesFile.write('\n')

#------------------------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------------------------

rawFile.close()
blastFile.close()
csvFile.close()
hitFile.close()
missFile.close()
failedFile.close()
uncountedFile.close()
speciesFile.close()
targetFile.close()
targetFileFa.close()


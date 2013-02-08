#!/usr/bin/python

# Copyright (c) 2011,2012,2013, Stephen Fisher and Junhyong Kim, University of
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
by: S. Fisher, 2011

usage: parseBlast.py <reads fasta file> <blast prefix>
assumes the input file
  <blast prefix>.txt
and will generate 
  <blast prefix>.csv - list of each hit per read and count of total number of hits per read
  <blast prefix>.hits - the alignments, when there is more than one alignment for a read
"""

import sys, os

#------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------

DEBUG = 0

if DEBUG: print 'DEBUG MODE: ON'

# expect 2 args
if len(sys.argv) < 3:
    print 'Usage: parseBlast.py <reads fasta file> <blast prefix>'
    print '  assumes the input file:'
    print '\t<blast prefix>.txt'
    print '  and will generate'
    print '\t<blast prefix>.csv - list of each hit per read and count of total number of hits per read'
    print '\t<blast prefix>.hits - the alignments, when there is more than one alignment for a read'
    print '\tfailed-<reads fasta file> - fasta file with reads that did not align'
    print '\n\tmouse.fa - fasta file with all mouse reads'
    print '\trat.fa - fasta file with all rat reads'
    print '\thuman.fa - fasta file with all human reads'
    print '\tbact.fa - fasta file with all bacteria reads'
    print '\tdanio.fa - fasta file with all zebrafish, danio rerio reads'
    print '\n\tmouse.fq - fastq file with all mouse reads'
    print '\trat.fq - fastq file with all rat reads'
    print '\thuman.fq - fastq file with all human reads'
    print '\tbact.fq - fastq file with all bacteria reads'
    print '\tdanio.fq - fastq file with all zebrafish, danio rerio reads'
    sys.exit()


RAW_FILE = sys.argv[1]
BLAST_FILE = sys.argv[2]
BLAST_PATH = os.path.dirname(RAW_FILE)
if len(BLAST_PATH) > 0:
    BLAST_PATH = BLAST_PATH + "/"

rawFile = open(RAW_FILE, 'r')
blastFile = open(BLAST_FILE+".txt", 'r')
csvFile = open(BLAST_FILE+".csv", 'w')
hitFile = open(BLAST_FILE+".hits", 'w')
failedFile = open(BLAST_PATH+"failed."+os.path.basename(RAW_FILE), 'w')

mouseFile = open(BLAST_PATH+'mouse.txt', 'w')
ratFile = open(BLAST_PATH+'rat.txt', 'w')
humanFile = open(BLAST_PATH+'human.txt', 'w')
bactFile = open(BLAST_PATH+'bact.txt', 'w')
danioFile = open(BLAST_PATH+'danio.txt', 'w')

mouseFileFa = open(BLAST_PATH+'mouse.fa', 'w')
ratFileFa = open(BLAST_PATH+'rat.fa', 'w')
humanFileFa = open(BLAST_PATH+'human.fa', 'w')
bactFileFa = open(BLAST_PATH+'bact.fa', 'w')
danioFileFa = open(BLAST_PATH+'danio.fa', 'w')

mouseFileFq = open(BLAST_PATH+'mouse.fq', 'w')
ratFileFq = open(BLAST_PATH+'rat.fq', 'w')
humanFileFq = open(BLAST_PATH+'human.fq', 'w')
bactFileFq = open(BLAST_PATH+'bact.fq', 'w')
danioFileFq = open(BLAST_PATH+'danio.fq', 'w')


def error(line):
    # this should never happen
    print "Got to this line in error: ", line, "\n"
    rawFile.close()
    blastFile.close()
    csvFile.close()
    hitFile.close()
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

def median(vals):
    if len(vals) == 0: return -1
    vals = sorted(vals)
    size = len(vals)
    if size % 2 == 1: return vals[(size - 1) / 2]
    else: return (vals[size/2 - 1] + vals[size/2]) / 2

numReads = 0
numFails = 0
numHits = 0

# we count each time a read is aligned to one of the following. Since
# we only want to count the species once per read, we use a flag to
# denote whether or not the read has already been counted.
numHuman = 0
numMouse = 0
numRat = 0
numBact = 0
numDanio = 0

numHsMm = 0
numHsRn = 0
numMmRn = 0
numHsMmRn = 0

# store values to compute median scores
humanTotId = []
humanPercId = []
humanPercGap = []
mouseTotId = []
mousePercId = []
mousePercGap = []
ratTotId = []
ratPercId = []
ratPercGap = []
bactTotId = []
bactPercId = []
bactPercGap = []
danioTotId = []
danioPercId = []
danioPercGap = []

while True:
    line = blastFile.readline()
    if not line: break

    if not 'Query=' in line:
        # keep going until fine a "Query="
        continue
    
    ###############################################
    # starting a new read
    query = line
    if DEBUG: print "\n---------------------------------------------------------------------------------\n", line,

    rawSeq, seqHeader = getRaw()
    #csvFile.write("Seq=   " + rawSeq)
    #csvFile.write(line)

    # count the read
    numReads += 1

    # skip next 5 lines
    for cnt in range(0, 5): line = getLine()

    # test if no alignments 
    if 'No hits' in line:
        if DEBUG: print "  Hits: 0"

        csvFile.write("0, " + rawSeq + ", " + seqHeader + "\n")

        # generate fasta file with sequences that didn't align
        failedFile.write(seqHeader + "\n" + rawSeq + "\n")
        
        # didn't align, so count the no hit and continue to next read
        numFails += 1

        # this will skip to next "Query=" statement
        continue
    
    ###############################################
    # process hits
    hits = 0
    humanFlag = 0
    mouseFlag = 0
    ratFlag = 0
    bactFlag = 0
    danioFlag = 0
    HsMmFlag = 0
    HsRnFlag = 0
    MmRnFlag = 0
    HsMmRnFlag = 0
    while True:
        line = getLine()

        # if empty line then ran out of alignments
        if line == "": break

        hits += 1

        lLine = line.lower() # test lowercase version of line
        if not humanFlag and (('human' in lLine) or ('sapiens' in lLine)): 
            numHuman += 1
            humanFlag = 1 # only count once per read
            humanFile.write(seqHeader + "\t" + rawSeq + "\t" + lLine + "\n")
            humanFileFa.write(seqHeader + "\n" + rawSeq + "\n")
            humanFileFq.write("@" + seqHeader + "\n" + rawSeq + "\n+\n" + "#"*101 + "\n")
            if mouseFlag and not HsMmFlag:
                numHsMm += 1
                HsMmFlag = 1
            if ratFlag and not HsRnFlag:
                numHsRn += 1
                HsRnFlag = 1
            if mouseFlag and ratFlag and not HsMmRnFlag:
                numHsMmRn += 1
                HsMmRnFlag = 1
        if not mouseFlag and (('mouse' in lLine) or ('musculus' in lLine)): 
            numMouse += 1
            mouseFlag = 1
            mouseFile.write(seqHeader + "\t" + rawSeq + "\t" + lLine + "\n")
            mouseFileFa.write(seqHeader + "\n" + rawSeq + "\n")
            mouseFileFq.write("@" + seqHeader + "\n" + rawSeq + "\n+\n" + "#"*101 + "\n")
            if humanFlag and not HsMmFlag:
                numHsMm += 1
                HsMmFlag = 1
            if ratFlag and not MmRnFlag:
                numMmRn += 1
                MmRnFlag = 1
            if humanFlag and ratFlag and not HsMmRnFlag:
                numHsMmRn += 1
                HsMmRnFlag = 1
        if not ratFlag and (('rattus' in lLine) or ('norvegicus' in lLine)): 
            numRat += 1
            ratFlag = 1
            ratFile.write(seqHeader + "\t" + rawSeq + "\t" + lLine + "\n")
            ratFileFa.write(seqHeader + "\n" + rawSeq + "\n")
            ratFileFq.write("@" + seqHeader + "\n" + rawSeq + "\n+\n" + "#"*101 + "\n")
            if mouseFlag and not MmRnFlag:
                numMmRn += 1
                MmRnFlag = 1
            if humanFlag and not HsRnFlag:
                numHsRn += 1
                HsRnFlag = 1
            if mouseFlag and humanFlag and not HsMmRnFlag:
                numHsMmRn += 1
                HsMmRnFlag = 1
        if not bactFlag and ('bacter' in lLine): 
            numBact += 1
            bactFlag = 1
            bactFile.write(seqHeader + "\t" + rawSeq + "\t" + lLine + "\n")
            bactFileFa.write(seqHeader + "\n" + rawSeq + "\n")
            bactFileFq.write("@" + seqHeader + "\n" + rawSeq + "\n+\n" + "#"*101 + "\n")
        if not danioFlag and (('zebrafish' in lLine) or ('danio rerio' in lLine)):
            numDanio +=1
            danioFlag = 1
            danioFile.write(seqHeader + "\t" + rawSeq + "\t" + lLine + "\n")
            danioFileFa.write(seqHeader + "\n" + rawSeq + "\n")
            danioFileFq.write("@" + seqHeader + "\n" + rawSeq + "\n+\n" + "#"*101 + "\n")

        # print alignment
        if DEBUG: print " ", line
        #csvFile.write("     " + line)

    if DEBUG: print "  Hits:", hits
    csvFile.write(str(hits) + ", " + rawSeq + ", " + seqHeader + "\n")

    # count number of hits
    if hits > 0: numHits += 1

    ###############################################
    # XXX If this is set to "<2" then only multiple alignments will be included in hits file
    # if no alignments, then start on next query
    if hits == 0: continue

    ###############################################
    # process multiple alignments

    # output alignment info to hitFile
    hitFile.write("------------------------------------------------------------------------------------------\n" )
    #hitFile.write("Seq=   " + rawSeq + "\n")
    hitFile.write(query)
    for hit in range(0, hits):
        #hitFile.write("\n  Alignment: " + str(hit+1) + "\n")

        # step through until we find an alignment
        while True:
            line = getLine()

            # stop when we get to the alignment
            if ">" in line: break

            # if we don't find any alignments, then there's a problem with the blast file
            if "Effective search space used" in line: error(line)

        # get alignement descriptor. It might be multiple lines and
        # will always begin with a ">" and end with a line with the length
        desc = ""
        while not "Length=" in line:
            desc += line + " "
            line = getLine()
        if len(desc) > 80: desc = desc[:77] + "..."
        hitFile.write("\n   " + desc + "\n")

        # skip blank line
        getLine()

        # print line with E value
        line = getLine()
        hitFile.write("  " + line + "\n")

        # print line with Identities and Gaps
        line = getLine()
        line = line.replace('(', '')
        line = line.replace('%)', ' ')
        line = line.replace('/', ' ')
        vals = line.split()
        desc= desc.lower() # test lowercase version of description
        if ('human' in desc) or ('sapiens' in desc): 
            humanTotId.append(int(vals[3]))
            humanPercId.append(int(vals[4]))
            humanPercGap.append(int(vals[10]))
        if ('mouse' in desc) or ('musculus' in desc): 
            mouseTotId.append(int(vals[3]))
            mousePercId.append(int(vals[4]))
            mousePercGap.append(int(vals[10]))
        if ('rattus' in desc): 
            ratTotId.append(int(vals[3]))
            ratPercId.append(int(vals[4]))
            ratPercGap.append(int(vals[10]))
        if ('bacter' in desc): 
            bactTotId.append(int(vals[3]))
            bactPercId.append(int(vals[4]))
            bactPercGap.append(int(vals[10]))
        if ('zebrafish' in desc) or ('danio rerio' in desc):
            danioTotId.append(int(vals[3]))
            danioPercId.append(int(vals[4]))
            danioPercGap.append(int(vals[10]))
                

        # skip next 2 lines
        for cnt in range(0, 2): getLine()

        # an alignment has 3 lines
        for cnt in range(0, 3): 
            line = getLine()
            hitFile.write("      " + line + "\n")
        
        # skip blank line
        line = getLine()

        # test if the alignment continues
        line = getLine()
        if "Query" in line:
            hitFile.write("\n")
            for cnt in range(0, 3): 
                hitFile.write("      " + line + "\n")
                line = getLine()

    hitFile.write("\n")


rawFile.close()
blastFile.close()
csvFile.close()
hitFile.close()
failedFile.close()

mouseFile.close()
humanFile.close()
ratFile.close()
bactFile.close()
danioFile.close()

mouseFileFa.close()
humanFileFa.close()
ratFileFa.close()
bactFileFa.close()
danioFileFa.close()

mouseFileFq.close()
humanFileFq.close()
ratFileFq.close()
bactFileFq.close()
danioFileFq.close()

print "\nNum reads:", numReads
print "Num reads that did not align:", numFails
print 'Number hits: %d ( %.1f%% )' % (numHits, (100*float(numHits)/float(numReads)))

print "\nRat Hits\t\t", numRat, "\tIdentity\t", median(ratTotId), "\tPerc Ident\t", median(ratPercId)
print "Mouse Hits\t\t", numMouse, "\tIdentity\t", median(mouseTotId), "\tPerc Ident\t", median(mousePercId)
print "Human Hits\t\t", numHuman, "\tIdentity\t", median(humanTotId), "\tPerc Ident\t", median(humanPercId)
print "Bacteria Hits\t\t", numBact, "\tIdentity\t", median(bactTotId), "\tPerc Ident\t", median(bactPercId)
print "Zebrafish Hits\t\t", numDanio, "\tIdentity\t", median(danioTotId), "\tPerc Ident\t", median(danioPercId)
print ""
print "Hits conserved Human and Mouse\t\t", numHsMm
print "Hits conserved Human and Rat\t\t", numHsRn
print "Hits conserved Mouse and Rat\t\t", numMmRn
print "Hits conserved Human, Mouse, and Rat\t", numHsMmRn,

#print "\nRat Hits\t\t", numRat, "\tIdentity\t", median(ratTotId), "\tPerc Ident\t", median(ratPercId), "\tPerc Gap\t", median(ratPercGap)
#print "Mouse Hits\t\t", numMouse, "\tIdentity\t", median(mouseTotId), "\tPerc Ident\t", median(mousePercId), "\tPerc Gap\t", median(mousePercGap)
#print "Human Hits\t\t", numHuman, "\tIdentity\t", median(humanTotId), "\tPerc Ident\t", median(humanPercId), "\tPerc Gap\t", median(humanPercGap)
#print "Bacteria Hits\t\t", numBact, "\tIdentity\t", median(bactTotId), "\tPerc Ident\t", median(bactPercId), "\tPerc Gap\t", median(bactPercGap)
#print "Zebrafish Hits\t\t", numDanio, "\tIdentity\t", median(danioTotId), "\tPerc Ident\t", median(danioPercId),"\tPerc Gap\t", median(danioPercGap)

if DEBUG:
    print humanTotId
    print humanPercId
    print humanPercGap
    print mouseTotId
    print mousePercId
    print mousePercGap
    print ratTotId
    print ratPercId
    print ratPercGap
    print bactTotId
    print bactPercId
    print bactPercGap
    print danioTotId
    print danioPercId
    print danioPercGap


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

usage: flyBase2RUM.py <exons GFF> <output>

LOCATION OF GFF FILE
ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.49_FB2013_01/gff/dmel-all-no-analysis-r5.50.gff.gz
PREPARING INPUT FILE:
  grep -P "exon\t" dmel-all-no-analysis-r5.50.gff > exons.tmp
  uniq exons.tmp > exons.gff
  rm exons.tmp

RUM requires that start is 0-based and end is 1-based. This is based
on a USCS requirement. The exons from FlyBase are 1-based, so the
start positions are adjusted as needed.

There are transcripts in the FlyBase transcriptome that contain exons
on both the + and - strands. The SINGLE_STRAND can be used to specify
how these transcripts should be handled.
  http://en.wikipedia.org/wiki/Trans-splicing

REQUIRES python 2.7: viewvalues()
"""

#------------------------------------------------------------------------------------------
# INITIALIZATIONS
#------------------------------------------------------------------------------------------

import sys

DEBUG = False
if DEBUG: print 'DEBUG MODE: ON'

# There are transcripts in the transcriptome that contain exons on
# both the + and - strands. When SINGLE_STRAND is true then the strand
# of a transcript will be the strand of the first exon for that
# transcript and any exons found on the other strand will be
# ignored. If SINGLE_STRAND is false then multi-stranded transcripts
# will be split into two transcripts (one for each strand) and the
# strings a ".1" and ".2" will be appended to the respective
# transcript names.
SINGLE_STRAND = False

# expect 2 argument
if len(sys.argv) < 3:
    print 'Usage: flyBase2RUM.py <exons GFF> <output>\n'
    print 'LOCATION OF GFF FILE'
    print 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.49_FB2013_01/gff/dmel-all-no-analysis-r5.50.gff.gz\n'
    print 'PREPARING INPUT FILE:'
    print '  grep -P "exon\\t" dmel-all-no-analysis-r5.50.gff > exons.tmp'
    print '  uniq exons.tmp > exons.gff'
    print '  rm exons.tmp\n'
    print 'RUM requires that start is 0-based and end is 1-based. This is based'
    print 'on a USCS requirement. The exons from FlyBase are 1-based, so the'
    print 'start positions are adjusted as needed.\n'
    print 'There are transcripts in the FlyBase transcriptome that contain exons'
    print 'on both the + and - strands. The SINGLE_STRAND can be used to specify'
    print 'how these transcripts should be handled.'
    sys.exit()

EXON_FILE = sys.argv[1]
OUT_FILE = sys.argv[2]

exonFile = open(EXON_FILE, 'r')
outFile = open(OUT_FILE, 'w')

#------------------------------------------------------------------------------------------
# DEFINE TRANSCRIPT CLASS
#------------------------------------------------------------------------------------------

class Transcript:
    # chrom is chromosome
    # strand is +/-
    # starts is 0-based string of exon start positions (RUM requirement based on UCSC requirement)
    # ends is 1-based string of exon end positions

    def __init__(self, name, chrom, strand, start, end):
        self.name, self.chrom, self.strand = name, chrom, strand
        # start is 1-based but needs to be 0-based
        start = int(start) - 1
        self.starts = str(start)
        # end is 1-based
        self.ends = end

        # if transcript contains exons on both +/- strands, then we
        # track both sets of exons. This comes into play when we find
        # an exon that is on the opposite strand from the first exon
        # in the transcript. A ".0" and ".1" will be appended to the
        # name of transcripts with two strands.
        self.starts0 = ""
        self.ends0 = ""
        

    def addExon(self, chrom, strand, start, end):
        # consistency checks for input
        if self.chrom != chrom:
            print "ERROR: transcript (" + self.name + ") spans two chromosomes."
            sys.exit()

        start = int(start) - 1
        if self.strand == strand:
            # add exon to transcript
            self.starts += "," + str(start)
            self.ends += "," + end
        else:
            # if transcript contains exons on both +/- strands, then we track both sets of exons.
            if SINGLE_STRAND:
                print "WARNING: Transcript (" + self.name + ") on two strands. Skipping exon - Chr: %s\tStrand: %s\tStart: %s\tEnd: %s" % (chrom, strand, start, end)
            else:
                print "WARNING: Transcript (" + self.name + ") on two strands. Splitting transcript - Chr: %s\tStrand: %s\tStart: %s\tEnd: %s" % (chrom, strand, start, end)
            if not SINGLE_STRAND:
                if self.starts0 == "":
                    self.starts0 = str(start)
                    self.ends0 = end
                else:
                    self.starts0 += "," + str(start)
                    self.ends0 += "," + end

    def writeToFile(self):
        # name   chrom   strand  exonStarts      exonEnds

        if SINGLE_STRAND or (self.starts0 == ""):
            # only one strand
            outFile.write("%s\t%s\t%s\t%s\t%s\n" % (self.name, self.chrom, self.strand, self.starts, self.ends))
        else:
            # print first strand
            outFile.write("%s\t%s\t%s\t%s\t%s\n" % (self.name + ".1", self.chrom, self.strand, self.starts, self.ends))

            # print second strand
            strand0 = "+" if (self.strand == "-") else "-"
            outFile.write("%s\t%s\t%s\t%s\t%s\n" % (self.name + ".2", self.chrom, strand0, self.starts0, self.ends0))
            

#------------------------------------------------------------------------------------------
# BUILD LIST OF EXONS
#------------------------------------------------------------------------------------------

transcripts = {}
for line in exonFile:
    # GFF file is tab delimited
    fields = line.split('\t')
    chrom = fields[0]
    start = fields[3]  # starting coordinate (1-based)
    end = fields[4]  # ending coordinate (1-based)
    strand = fields[6] # feature strand
    attribs = fields[8].split(';') # semicolon separated list of attributes (eg ID=FBgn0031208:2;Parent=FBtr0300690)
    
    # get transcript ID(s) from attributes. Each exon might be in multiple transcripts
    # ex: ID=FBgn0003963:13;Name=ush:13;Parent=FBtr0078063,FBtr0329895;parent_type=mRNA
    for attrib in attribs:
        # the Parent key is the transcript ID. Skip other attributes
        if 'Parent=' not in attrib: continue

        # get list of transcript IDs
        # ex: Parent=FBtr0078063,FBtr0329895
        ids = (attrib.split('='))[1].split(',')
    
        for name in ids:
            # add exon to transcript dictionary
            if name in transcripts:
                transcripts[name].addExon(chrom, strand, start, end)
            else:
                transcripts[name] = Transcript(name, chrom, strand, start, end)

#------------------------------------------------------------------------------------------
# GENERATE OUTPUT FILE
#------------------------------------------------------------------------------------------

# this header is needed by RUM
outFile.write("#name\tchrom\tstrand\texonStarts\texonEnds\n")

cnt = 0
# viewvalues() requires python 2.7
for transcript in transcripts.viewvalues():
    transcript.writeToFile()

#------------------------------------------------------------------------------------------
# CLEAN UP
#------------------------------------------------------------------------------------------

exonFile.close()
outFile.close()

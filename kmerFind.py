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
by: S. Fisher, J. Shallcross
usage: kmerTrim.py [-h] [-v] [-c] -f FIRST_FQ [-r SECOND_FQ] -o OUTPUT_PREFIX

Return all reads that contain the specified feature, allowing for one mismatch. The k-mer is removed from the read if the -c flag is present.
"""

#xTODO: create (or access) directories for each barcode
#xTODO: pass list of barcode regex to multithreaded portion
#xTODO: change kmer from global, then call mapReads for each kmer
#xTODO: find a way to mark which kmer mapped, accounting for multiple maps and multiple cutting points 
#xTODO: process reads, work out if it maps to one, many, or no kmers, output accordingly 
#xTODO: remove/minimize globals as much as possible

#------------------------------------------------------------------------------------------
# INITIALIZATIONS
#------------------------------------------------------------------------------------------

import sys, os, errno, argparse, re, multiprocessing
#Partial w/ multiprocessing requires python 2.7+. Used to reduce global variable load and b/c nested functions don't work w/ multiprocessing since they can't be pickled
#If nessesary it should still work in 2.6 if single threaded (map rather than pool.map)
from functools import partial 
from Bio.Seq import Seq

import matplotlib
matplotlib.use('Agg')  #for headless use
import matplotlib.pyplot as plt

DEBUG = True
if DEBUG: sys.stderr.write('DEBUG MODE: ON\n')

VERSION = '0.6.2'

# indecies for the read set
HEADER = 'header'
SEQUENCE = 'seq'
QUALS = 'quals'
MAP_LOCATION = 'mapLoc'
MAPPED = 'mapped'
DIRECTION = 'direction'
LOCATION = 'location'
FOUND_SEQ = 'foundSeq'

CHUNK_SIZE = 5000

# print error and quit
def quitOnError(msg):
    print 'ERROR:', msg
    sys.exit(0)

#------------------------------------------------------------------------------------------
# SET UP SEARCH TERMS
#------------------------------------------------------------------------------------------

# k-mer without any mismatches
#KMER = "" #"GTAGAGTTTTTT"
#KMER_RC = "" #"AAAAAACTCTAC"

# k-mer allowing for 1 missmatch
#KMER_1_MISMATCH = re.compile('([ATCG]TAGAGTTTTTT)|(G[ATCG]AGAGTTTTTT)|(GT[ATCG]GAGTTTTTT)|(GTA[ATCG]AGTTTTTT)|(GTAG[ATCG]GTTTTTT)|(GTAGA[ATCG]TTTTTT)|(GTAGAG[ATCG]TTTTT)|(GTAGAGT[ATCG]TTTT)|(GTAGAGTT[ATCG]TTT)|(GTAGAGTTT[ATCG]TT)|(GTAGAGTTTT[ATCG]T)|(GTAGAGTTTTT[ATCG])')

#KMER_1_MISMATCH_RC = re.compile('([ATCG]AAAAACTCTAC)|(A[ATCG]AAAACTCTAC)|(AA[ATCG]AAACTCTAC)|(AAA[ATCG]AACTCTAC)|(AAAA[ATCG]ACTCTAC)|(AAAAA[ATCG]CTCTAC)|(AAAAAA[ATCG]TCTAC)|(AAAAAAC[ATCG]CTAC)|(AAAAAACT[ATCG]TAC)|(AAAAAACTC[ATCG]AC)|(AAAAAACTCT[ATCG]C)|(AAAAAACTCTA[ATCG])')


def kmerRegEx(seq):
    regex = ""
    for i in xrange(len(seq)):
	if i != 0: 
	    regex += '|'
	regex += '(' + seq[:i] + '[ATCG]' + seq[i+1:] + ')'
    #regex = re.compile(regex)
    return regex

def idxToKmer(index):
    illuminaIndexes = { 
	    1 : "ATCACG",
	    2 : "CGATGT",
	    3 : "TTAGGC",
	    4 : "TGACCA",
	    5 : "ACAGTG",
	    6 : "GCCAAT",
	    7 : "CAGATC",
	    8 : "ACTTGA",
	    9 : "GATCAG",
	    10 : "TAGCTT",
	    11 : "GGCTAC",
	    12 : "CTTGTA",
	    13 : "AGTCAA",
	    14 : "AGTTCC",
	    15 : "ATGTCA",
	    16 : "CCGTCC",
	    17 : "GTAGAG",
	    18 : "GTCCGC",
	    19 : "GTGAAA",
	    20 : "GTGGCC",
	    21 : "GTTTCG",
	    22 : "CGTACG",
	    23 : "GAGTGG",
	    24 : "GGTAGC",
	    25 : "ACTGAT",
	    26 : "ATGAGC",
	    27 : "ATTCCT", } 
    if index not in illuminaIndexes:
	raise ValueError("%d not a valid Illumina barcode index" % index)
    return illuminaIndexes[index]

def kmerToIdx(kmer):
    illuminaIndexes = { 
	   "ATCACG" : 1,
	   "CGATGT" : 2,
	   "TTAGGC" : 3,
	   "TGACCA" : 4,
	   "ACAGTG" : 5,
	   "GCCAAT" : 6,
	   "CAGATC" : 7,
	   "ACTTGA" : 8,
	   "GATCAG" : 9,
	   "TAGCTT" : 10,
	   "GGCTAC" : 11,
	   "CTTGTA" : 12,
	   "AGTCAA" : 13,
	   "AGTTCC" : 14,
	   "ATGTCA" : 15,
	   "CCGTCC" : 16,
	   "GTAGAG" : 17,
	   "GTCCGC" : 18,
	   "GTGAAA" : 19,
	   "GTGGCC" : 20,
	   "GTTTCG" : 21,
	   "CGTACG" : 22,
	   "GAGTGG" : 23,
	   "GGTAGC" : 24,
	   "ACTGAT" : 25,
	   "ATGAGC" : 26,
	   "ATTCCT" : 27, } 
    if kmer not in illuminaIndexes:
	raise ValueError("%s not a valid Illumina barcode sequence" % kmer)
    return illuminaIndexes[kmer]

# number of bases to shift index in order to mark 3' end of k-mer, rather than 5' end of k-mer
#KMER_OFFSET = 0

# only find k-mers that are within specified distance from end of
# read. The values below are relative to the 5' end of the read and
# the 5' end of the k_mer. The values below assume the k_mer will be
# at either the 3' end of the first read or the 5' end of the second
# read.
#THRESH_1 = 0 #45 # for first read k_mer 5' location > THRESH_1
#THRESH_2 = 0 #10 # for second read k_mer 5' location < THRESH_2

#------------------------------------------------------------------------------------------
# FUNCTIONS FOR HANDLING READING AND WRITING OF READS
#------------------------------------------------------------------------------------------

def nextRead(inFile):
    """
    Returns a set consisting of each line from the read, a length
    indicator and a mapped flag. Returns empty set if no more reads.
    """
    line = inFile.readline().rstrip()
    if len(line) == 0:
        return {}

    read = {}
    try:
        # load read values from fastq file
        read[HEADER] = line
        # make sure sequence is upper case for comparisons
        read[SEQUENCE] = inFile.readline().rstrip().upper()
        # this is the + which we can ignore for now

        # XXX CONFIRM THAT THE + IN THE FASTQ FILE IS A CONSTANT.
        inFile.readline()

        read[QUALS] = inFile.readline().rstrip()
        
        # store additional read information in read set
        #read[MAP_LOCATION] = -1
        #read[MAPPED] = False

        return read
    except:
        msg = 'Problem loading read:', line
        quitOnError(msg)

def writeRead(read, outFile):
    """
    Writes a read to the read file. Read is expected to be a 4-item
    list as returned by nextRead().
    """
    outFile.write(read[HEADER] + '\n') # header
    outFile.write(read[SEQUENCE] + '\n') # sequence
    outFile.write('+\n') # +
    outFile.write(read[QUALS] + '\n') # quals
    
#------------------------------------------------------------------------------------------
# MAPPING FUNCTION
#------------------------------------------------------------------------------------------

def kmerFound(index, isFirstRead, t1, t2):
    """
    Return True if kmer found within proper position on read.
    """

    # test this first to speed up processing
    if index == -1: 
	return False

    if isFirstRead:
        if index > t1 or t1 == 0: return True
    else:
        if index < t2 or t2 == 0: return True
        #if (index > -1) and (index < THRESH_2): return True

    # if we get this far then either we found a match but in the wrong location.
    return False


def mapReads(read, isFirstRead, bcExp, cut_reads=False):
    """
    Map k-mer to a read using regex, cutting the read as required.
    """
    seq = read[SEQUENCE]
    quals = read[QUALS]
    kmer = bcExp["KMER"]
    kmer_rc = bcExp["KMER_RC"]
    t1 = bcExp["T1"]
    t2 = bcExp["T2"]

    read[kmer] = {}
    read[kmer][MAPPED] = False
    read[kmer][MAP_LOCATION] = -1
    read[kmer][FOUND_SEQ] = '-'

    # first look for full string. that will be faster than regex
    index = seq.find(kmer)

    if kmerFound(index, isFirstRead, t1, t2):
	read[kmer][DIRECTION] = "forward"
	read[kmer][FOUND_SEQ] = kmer
    else: # didn't find k-mer so check reverse compliment
        index = seq.find(kmer_rc)

	if kmerFound(index, isFirstRead, t1, t2):
	    read[kmer][DIRECTION] = "reverse"
	    read[kmer][FOUND_SEQ] = kmer_rc
	else: 
	    # still didn't find k-mer so check regular expression
	    regexOut = bcExp["KMER_REGEX"].search(seq)
	    if regexOut:
		index = regexOut.start()

	    if kmerFound(index, isFirstRead, t1, t2):
		read[kmer][DIRECTION] = "forward"
		read[kmer][FOUND_SEQ] = regexOut.group(0) #returns the matching string right?
	    else:
		# still didn't find k-mer so check regular expression reverse compliment
		regexOut = bcExp["KMER_RC_REGEX"].search(seq)
		if regexOut:
		    index = regexOut.start()

		if kmerFound(index, isFirstRead, t1, t2):
		    read[kmer][DIRECTION] = "reverse"
		    read[kmer][FOUND_SEQ] = regexOut.group(0) #returns the matching string right?

    if kmerFound(index, isFirstRead, t1, t2):
        # we found the string so trim read if required
        read[kmer][MAPPED] = True
        read[kmer][MAP_LOCATION] = index
        
        if cut_reads == "cut":
            if isFirstRead:
                # mate 1 is trimmed from index location to 3' end
                cutPos = index

                # leave at least 1 base
                if cutPos == 0: cutPos = 1
                    
                # trim read
                seq = seq[:cutPos]
                quals = quals[:cutPos]
            else:
                # 5' trimming. shift index to the 3' end of the k-mer
                cutPos = index + len(kmer)

                # leave at least 1 base
                if cutPos == len(seq): cutPos = cutPos - 1

                # trim read
                seq = seq[cutPos:]
                quals = quals[cutPos:]


	    # if in debug mode then output the region of the read that we are cutting off
	    if DEBUG: read[kmer]["cutPos"] = seq[cutPos:] + '\n'
	
	# Instead of cutting reads, just mask bc sequence with Ns
	elif cut_reads == "mask":
	    # ie 123456789 => 123NNN789 if kmer=456
	    seq = seq[:index] + 'N' * len(kmer) + seq[index+len(kmer):]
	    quals = quals[:index] + '#' * len(kmer) + quals[index+len(kmer):]

    read[kmer][SEQUENCE] = seq
    read[kmer][QUALS] = quals

    return read
    
def mapHelper(pair, barcodes, cut_reads):
    for bcSet in barcodes:
	firstRead = pair[0]
	secondRead  = pair[1]
	firstRead = mapReads(firstRead, True, bcSet, cut_reads)
	if secondRead: secondRead = mapReads(secondRead, False, bcSet, cut_reads)

    return (firstRead, secondRead)

#------------------------------------------------------------------------------------------
# OPEN READ INPUT AND OUTPUT FILES
#------------------------------------------------------------------------------------------
def main():


    argParser = argparse.ArgumentParser(version=VERSION, 
					description='Find k-mer in NGS reads.',
					formatter_class=argparse.RawDescriptionHelpFormatter,
					epilog='' +
					'Return all reads that contain the specified feature, allowing for one mismatch. The k-mer is removed from the read.\n'
    )


    argParser.add_argument( '-b', dest='barcode', action='store', required=True,
			    help='Barcode sequence(s) or Illumina index(es) to be selected for. Multiple barcodes should be comma (,) separated.' )
    argParser.add_argument( '-customExp', dest='isRegex', action='store_true', default=False,
			    help='Treat [barcode] as a custom python regex. Must be a valid regex as of 2.7.6, and must contain at least one capturing group (default: no)' )
    argParser.add_argument( '-rx2', dest='regex2', action='store', default='X^', 
			    help='Second regex used eith customExp as the reverse read. Defaults to no matching.' )
    argParser.add_argument( '-c', dest='cutReads', action='store_true', default=False,
			    help='Reads containing the k-mer will be trimmed as follows: (1) for mate-1 the read will be trimmed from the 5\' location of the k-mer to the 3\' end of the read, and (2) for mate-2 the read will be trimmed from the 3\' end of the k-mer to the 5\' end of the read. (default: no)' )
    argParser.add_argument( '-m', dest='maskReads', action='store_true', default=False,
			    help='Reads containing the k-mer will have the k-mer sequence masked with N\'s, but maintain the sequence on both sides of the k-mer.' )
    argParser.add_argument( '-f', dest='first_fq', action='store', required=True,
			    help='fastq file with reads to be mapped' )
    argParser.add_argument( '-r', dest='second_fq', 
			    help='second read with paired-end reads. This file is not present for single-end reads.' )
    argParser.add_argument( '-TF', dest='forward_thresh', type=int, default=0,
			    help="for first read k_mer 5' location > TF" )
    argParser.add_argument( '-TR', dest='reverse_thresh', type=int, default=0,
			    help="for second read k_mer 5' location > TR" )
    argParser.add_argument( '-t', dest='threads', type=int, default=1,
			    help='number of mapping threads.' )
    argParser.add_argument( '-prefix', dest='prefix', action='store', required=False, default='',
			    help='Prefix sequence, to be added to the front of all barcodes.' )
    argParser.add_argument( '-suffix', dest='suffix', action='store', required=False, default='',
			    help='Suffix sequence, to be added to the end of all barcodes.' )
    argParser.add_argument( '-psMode', dest='psMode', choices=['and', 'or'], action='store', required=False, default='and',
			    help='How to use the prefix/suffix. If "and", matched sequence must have both. If "or", if may have either the prefix or the suffix, or both.' ) 
    argParser.add_argument( '-o', dest='output_prefix', action='store', required=True,
			    help='prefix for output file(s). A \'_1\' will be appended to the first reads output file and if paired reads then a \'_2\' will be appended to the second reads file. The header for reads will be adjusted to include "::L:X" where X is the 5\' position of the k-mer within the read. If the -c flag is present then reads containing the k-mer will be trimmed as follows: (1) for mate-1 the read will be trimmed from the 5\' location of the k-mer to the 3\' end of the read and (2) for mate-2 the read will be trimmed from the 3\' end of the k-mer to the 5\' end of the read.' )

    clArgs = argParser.parse_args()
    if DEBUG: print clArgs

    # flag if paired-end reads
    PAIRED = False
    if clArgs.second_fq:
	PAIRED = True

    # flag if trimming reads
    CUT_READS = False
    if clArgs.cutReads:
	CUT_READS = "cut"
    elif clArgs.maskReads:
	CUT_READS = "mask"

    first_fq = clArgs.first_fq
    output_prefix = clArgs.output_prefix
    second_fq = clArgs.second_fq
    nThreads = clArgs.threads

    t1 = clArgs.forward_thresh 
    t2 = clArgs.reverse_thresh

    #body(first_fq, second_fq, output_prefix, threads, barcodes)
    #def body(first_fq, second_fq, output_prefix, nThreads, barcode):

    bcExps = []

    if not clArgs.isRegex: #is normal barcode
	barcodes = map(lambda b: idxToKmer(int(b)) if b.isdigit() else b, clArgs.barcode.split(',')) #Throws customized ValueError, but just let it pass through
	for b in barcodes:
	    #TODO allow for numeric BCs
	    try:
		bcIdx = kmerToIdx(b)

	    except ValueError:
		print "note: %s is not a known Illumina barcode" % b
		bcIdx = b

	    if clArgs.psMode == "and":
		b = clArgs.prefix + b + clArgs.suffix
		bRC = str(Seq(b).reverse_complement())
		bcExps.append({
		    "ILLUMINA_INDEX" : bcIdx, 
		    "KMER" : b, 
		    "KMER_RC" : bRC, 
		    "KMER_REGEX" : re.compile(kmerRegEx(b)),
		    "KMER_RC_REGEX" : re.compile(kmerRegEx(bRC)), 
		    "T1" : t1,
		    "T2" : t2,
		}) # Including thresholds t1/t2 in the bc expressions makes it easy if we ever want different thresholds for each bc
	    elif clArgs.psMode == "or":
		p_b = clArgs.prefix + b
		b_s = b + clArgs.suffix
		bcRegex = re.compile(kmerRegEx(p_b) + '|' + kmerRegEx(b_s)) # Combine to form a regex that will hit to either prefix + b or b + suffix
		bcName = p_b + ' | ' + b_s 
		# NOTE: "bcName" negates the gains of using string.find for this string in mapReads, since it should never be found. 
		# The current form, however, makes printing/tracking results easier. If "or" becomes the standard form, the string.find calls should be eliminated
		# On the other hand, the current form (0.5.2) is I/O bound, not CPU bound, when using multiple threads, so it may not matter
		# The other option is to use prefix + b + suffix, as above, which would find the rare case with all three parts exactly. 

		p_bRC = str(Seq(p_b).reverse_complement())
		b_sRC = str(Seq(b_s).reverse_complement())
		bcRegexRC = re.compile(kmerRegEx(p_bRC) + '|' + kmerRegEx(b_sRC))
		bcNameRC = p_bRC + '|' + b_sRC
		bcExps.append({
		    "ILLUMINA_INDEX" : bcIdx, 
		    "KMER" : bcName, 
		    "KMER_RC" : bcNameRC, 
		    "KMER_REGEX" : bcRegex,
		    "KMER_RC_REGEX" : bcRegexRC, 
		    "T1" : t1,
		    "T2" : t2,
		})
    else: #is custom regex
	barcodes = []


    # open input file
    firstReadIn = ''
    try: 
	firstReadIn = open(first_fq, 'r')
    except: 
	msg = 'Unable to load reads file ' + first_fq
	quitOnError(msg)

    # open unmapped file
    
    if len(barcodes) > 1:
	output_path = output_prefix + "_Unmapped/barcode.trim/"
    else:
	output_path = output_prefix

    try:
	os.makedirs(output_path)
    except OSError as exc:  # Python >2.5
	if exc.errno == errno.EEXIST and os.path.isdir(output_path): #Directory already exists
	    pass
	else:	#some other error
	    raise 


    try: 
	firstUnmapped = open(output_path + "unmapped_1.fq", 'w')
    except: 
	msg = 'Unable to open output file ' + output_path + "unmapped_1.fq"
	quitOnError(msg)


    if PAIRED:
	secondReadIn = ''
	try: 
	    secondReadIn = open(second_fq, 'r')
	except: 
	    msg = 'Unable to load reads file ' + second_fq
	    quitOnError(msg)
	try: 
	    secondUnmapped = open(output_path + "unmapped_2.fq", 'w')
	except: 
	    msg = 'Unable to open output file ' + output_prefix + output_midfix + "unmapped_2.fq"
	    quitOnError(msg)

    nTotalReadPairs = 0 # total number of read pairs processed
    nAmbiguous = 0
    stats = {}
    files = {}
    for bx in bcExps:
	b = bx["KMER"]
	bcIdx = bx["ILLUMINA_INDEX"]

	stats[b] = {}
	stats[b]["idx"] = bcIdx
	# track mapping stats
	stats[b]["nBothMapped"] = 0 # total number of read pairs in which both reads were mapped (equal to nFirstMapped if single-end)
	stats[b]["nFirstMapped"] = 0 # num first reads mapped
	stats[b]["nSecondMapped"] = 0 # num second reads mapped
	stats[b]["firstForwardMappedPos"] = []
	stats[b]["firstReverseMappedPos"] = []
	stats[b]["secondForwardMappedPos"] = []
	stats[b]["secondReverseMappedPos"] = []

	
	if len(barcodes) > 1:
	    output_path = output_prefix + '_' + str(bcIdx) + "/barcode.trim/"
	else:
	    output_path = output_prefix
	
	try:
	    os.makedirs(output_path)
	except OSError as exc:  # Python >2.5
	    if exc.errno == errno.EEXIST and os.path.isdir(output_path): #Directory already exists
		pass
	    else:	#some other error
		raise 

	# open file that will store non-discarded reads
	files[b] = {}
	files[b]["firstReadOut"] = ''
	try: 
	    files[b]["firstReadOut"] = open(output_path + "unaligned_1.fq", 'w')
	except: 
	    msg = 'Unable to open output file ' + output_path + "unaligned_1.fq"
	    quitOnError(msg)

	if PAIRED:

	    files[b]["secondReadOut"] = ''
	    try: 
		files[b]["secondReadOut"] = open(output_path + "unaligned_2.fq", 'w')
	    except: 
		msg = 'Unable to open output file ' + output_path + "unaligned_2.fq"
		quitOnError(msg)

	if DEBUG:
	    files[b]["firstReadDebug"] = ''
	    try: 
		files[b]["firstReadDebug"] = open(output_path + "unaligned_1.cut", 'w')
	    except: 
		msg = 'Unable to open output file ' + output_path + "unaligned_1.cut"
		quitOnError(msg)
	    if PAIRED:
		files[b]["secondReadDebug"] = ''
		try: 
		    files[b]["secondReadDebug"] = open(output_path + "unaligned_2.cut", 'w')
		except: 
		    msg = 'Unable to open output file ' + output_path + "unaligned_2.cut"
		    quitOnError(msg)

    #------------------------------------------------------------------------------------------
    # MAP READS.  LOAD ONE READ FROM BOTH FQ FILES AT SAME TIME. PROCESS
    # READ PAIR TOGETHER AND THEN EITHER WRITE TO OUTPUT FILE OR DISCARD.
    # ------------------------------------------------------------------------------------------


    xMapHelper = partial(mapHelper, barcodes=bcExps, cut_reads=CUT_READS) # partial allows for passing of fixed arguments, only chunk changes
    pool = multiprocessing.Pool(nThreads)

    #--------------------------------------------------------------------------------------
    # main loop
    done = False
    runningJob = None
    results = []

    while runningJob or not done:
	
	#--------------------------------------------------------------------------------------
	# build read chunk
	chunk = []
	if not done: #done set to True when we read the last record
	    for i in xrange(CHUNK_SIZE * nThreads):
		firstRead = nextRead(firstReadIn)
		if PAIRED: 
		    secondRead = nextRead(secondReadIn)
		else:
		    secondRead = None

		if firstRead == {}:
		    firstReadIn.close()
		    #firstReadOut.close()
		    if PAIRED: 
			secondReadIn.close()
		    done = True
		    break
			
		nTotalReadPairs += 1

		readPair = (firstRead, secondRead)
		chunk.append(readPair)

	#--------------------------------------------------------------------------------------
	# do mapping

	if nThreads == 1:
	    # This removes the overhead of multiprocessing while keeping the same structure
	    # Very useful for debugging, since stack traces don't get passed upwards from child processes
	    results = map(xMapHelper,chunk)
	else:
	    if runningJob:
		results = runningJob.get() #grab results, waiting if they're not done
		runningJob = None

	    if chunk: # reads left to map
		# map_async allows us to do stats and IO work while the subprocesses run
		runningJob = pool.map_async(xMapHelper, chunk) #start new running job

	    if done and not runningJob:
	        pool.close()
	        pool.join()

	#--------------------------------------------------------------------------------------
	# process output
	for firstRead, secondRead in results: #First iteration results is [] so this skips
	    
	    firstResults = set()
	    secondResults = set()

	    for bx in bcExps:
		b = bx["KMER"]
		if firstRead[b][MAPPED]:
		    firstResults.add(b)
	    if PAIRED:
		for bx in bcExps:
		    b = bx["KMER"]
		    if secondRead[b][MAPPED]:
			secondResults.add(b)

	    #--------------------------------------------------------------------------------------
	    # loop if neither read mapped to any barcode
	    if PAIRED: 
		if not firstResults and not secondResults: #no results for either
		    writeRead(firstRead, firstUnmapped)
		    writeRead(secondRead, secondUnmapped)
		    continue
	    elif not firstResults:
		writeRead(firstRead, firstUnmapped)
		continue

	    
	    #--------------------------------------------------------------------------------------
	    # check for ambiguous reads 
	    mappedBCs = firstResults.union(secondResults)
	    if len(mappedBCs) > 1: 
		nAmbiguous += 1
		#TODO add mapped kmers to header
		firstRead[HEADER] += "\tambig: " + str(list(mappedBCs))
		writeRead(firstRead, firstUnmapped)
		if PAIRED:
		    secondRead[HEADER] += "\tambig: " + repr(list(mappedBCs))
		    writeRead(secondRead, secondUnmapped)
		continue


	    #--------------------------------------------------------------------------------------
	    # compute trimming stats
	    # if we're here, the read mapped once & only once
	    for bx in bcExps:
		b = bx["KMER"]
		mapFirst = False
		mapSecond = False
		if firstRead[b][MAPPED]: 
		    stats[b]["nFirstMapped"] += 1
		    mapFirst = True
		    if firstRead[b][DIRECTION] == "forward":
			stats[b]["firstForwardMappedPos"].append(firstRead[b][MAP_LOCATION])
		    if firstRead[b][DIRECTION] == "reverse":
			stats[b]["firstReverseMappedPos"].append(firstRead[b][MAP_LOCATION])
		if PAIRED: 
		    if secondRead[b][MAPPED]: 
			stats[b]["nSecondMapped"] += 1
			mapSecond = True
			if secondRead[b][DIRECTION] == "forward":
			    stats[b]["secondForwardMappedPos"].append(secondRead[b][MAP_LOCATION])
			if secondRead[b][DIRECTION] == "reverse":
			    stats[b]["secondReverseMappedPos"].append(secondRead[b][MAP_LOCATION])
		    if mapFirst and mapSecond:
			# count mapped pair
			stats[b]["nBothMapped"] += 1
		elif mapFirst:
		    # not paired, so count single read as pair
		    stats[b]["nBothMapped"] += 1

		if DEBUG and "cutPos" in firstRead[b]: 
		    files[b]["firstReadDebug"].write(firstRead[b]["cutPos"])
		if PAIRED and DEBUG and "cutPos" in secondRead[b]:
		    files[b]["secondReadDebug"].write(secondRead[b]["cutPos"])
		#--------------------------------------------------------------------------------------
		# add mapping location to read header
		kmer_found = "[%s, %s]" % (firstRead[b][FOUND_SEQ], secondRead[b][FOUND_SEQ])
		firstRead[HEADER] += '::L:%d\tBC:%s' % (firstRead[b][MAP_LOCATION], kmer_found)
		if PAIRED: secondRead[HEADER] += '::L:%d\tBC:%s' % (secondRead[b][MAP_LOCATION], kmer_found)


		# CUT_READS can be False, "cut", or "mask". If cut or mask, the base sequence must be replaced with the altered 
		# sequence for that kmer
		if CUT_READS:
		    firstRead[SEQUENCE] = firstRead[b][SEQUENCE]
		    firstRead[QUALS] = firstRead[b][QUALS]
		    if PAIRED:
			secondRead[SEQUENCE] = secondRead[b][SEQUENCE]
			secondRead[QUALS] = secondRead[b][QUALS]


		#--------------------------------------------------------------------------------------
		# write read(s) to output file if mapped to this bc
		if mapFirst or mapSecond:
		    writeRead(firstRead, files[b]["firstReadOut"])
		    if PAIRED: writeRead(secondRead, files[b]["secondReadOut"])

    #End while
    for b in files:
	for f in files[b].values():
	    f.close()

    #firstReadOut.close()
    #secondReadOut.close()

    #------------------------------------------------------------------------------------------
    # OUTPUT TRIMMING STATS
    #------------------------------------------------------------------------------------------
    for bc, s in stats.items():
	b = str(s["idx"])
	if len(stats) > 1:
	    try: 
		statsfile = open(output_prefix + '_' + b + "/barcode.trim/stats.txt", 'w')
	    except: 
		msg = 'Unable to open output file ' + output_prefix + '_' + b + "/barcode.trim/stats.txt"
		quitOnError(msg)
	else:
	    statsfile = sys.stdout

	statsfile.write('Read pairs processed: %s\n' % nTotalReadPairs )
	statsfile.write('Kmer sequence: %s\n' % bc )
	statsfile.write('\nBoth first and second reads mapped: %s\n' % s["nBothMapped"])
	statsfile.write('\tFirst reads mapped: %s\n' % s["nFirstMapped"])
	statsfile.write('\tSecond reads mapped: %s\n\n' % s["nSecondMapped"])

	avgFirst = sum(s["firstForwardMappedPos"] + s["firstReverseMappedPos"]) / float(len(s["firstForwardMappedPos"] + s["firstReverseMappedPos"]))
	avgSecond = sum(s["secondForwardMappedPos"] + s["secondReverseMappedPos"]) / float(len(s["secondForwardMappedPos"] + s["secondReverseMappedPos"]))

	# tab delimited output to facilitate adding stats to compilation file
	fields = '\nnTotalReadsBeforeMapping\ttotalBothMapped\ttotalFirstMapped\ttotalSecondMapped\tfirstKmer\tfirstRC\tsecondKmer\tsecondRC\tavgPosFirst\tavgPosSecond\n'
	counts = '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n' % (nTotalReadPairs, s["nBothMapped"], s["nFirstMapped"], s["nSecondMapped"], len(s["firstForwardMappedPos"]), len(s["firstReverseMappedPos"]), len(s["secondForwardMappedPos"]), len(s["secondReverseMappedPos"]), avgFirst, avgSecond)

	statsfile.write(fields)
	statsfile.write(counts)

	statsfile.close()

	# heuristic attempt to get resonable bins
	# if left to it's own devices pyplot does weird things like putting a small bin in the middle of two larger bins
	# this should always put the small bin at the end, if one is needed. 
	# this is still awkward if say 49 bases were covered
	pMax = max(s["firstForwardMappedPos"] + s["firstReverseMappedPos"]) 
	pMin = min(s["firstForwardMappedPos"] + s["firstReverseMappedPos"])
	nPos = pMax - pMin + 1
	if nPos > 20 and nPos % 2 == 0: 
	    bins = range(pMin, pMax, 2)
	elif (nPos > 30 and nPos % 3 == 0) or nPos > 50:
	    bins = range(pMin, pMax, 3)
	else:
	    bins = range(pMin, pMax, 1)

	if len(stats) > 1:
	    figpath = output_prefix + '_' + b + "/barcode.trim/" + b
	else:
	    figpath = output_prefix + "unaligned"

	if len(s["firstForwardMappedPos"]) > 0: 
	    plt.hist(s["firstForwardMappedPos"], bins, alpha=0.5, label='kmer')
	if len(s["firstReverseMappedPos"]) > 0:
	    plt.hist(s["firstReverseMappedPos"], bins, alpha=0.5, label='kmerRC')
	if len(s["firstForwardMappedPos"]) > 0 or len(s["firstReverseMappedPos"]) > 0:
	    plt.legend(loc='upper right')
	    plt.title("First Read (%s)" % bc )
	    plt.savefig(figpath + "_1.hist.pdf")
	    plt.close()

	if PAIRED and len(s["secondForwardMappedPos"]) > 0 or len(s["secondReverseMappedPos"]) > 0:
	    pMax = max(s["secondForwardMappedPos"] + s["secondReverseMappedPos"]) 
	    pMin = min(s["secondForwardMappedPos"] + s["secondReverseMappedPos"])
	    nPos = pMax - pMin + 1
	    if nPos > 20 and nPos % 2 == 0: 
		bins = range(pMin, pMax, 2)
	    elif (nPos > 30 and  nPos % 3 == 0) or nPos > 50:
		bins = range(pMin, pMax, 3)
	    else:
		bins = range(pMin, pMax, 1)

	    if len(s["secondForwardMappedPos"]) > 0:
		plt.hist(s["secondForwardMappedPos"], bins, alpha=0.5, label='kmer')
	    if len(s["secondReverseMappedPos"]) > 0:
		plt.hist(s["secondReverseMappedPos"], bins, alpha=0.5, label='kmerRC')
	    plt.legend(loc='upper right')
	    plt.title("Second Read (%s)" % bc)
	    plt.savefig(figpath + "_2.hist.pdf")
	    plt.close()

if __name__ == "__main__":
    main()

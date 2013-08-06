#!/bin/bash

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

##########################################################################################
# SINGLE-END READS:
# INPUT: $SAMPLE/raw/unaligned_1.fq
# OUTPUT: $SAMPLE/trimAT/unaligned_1.fq, $SAMPLE/trimAdapters.stats.txt, $SAMPLE/trimPolyAT.stats.txt
#         intermediate file $SAMPLE/trimAD/unaligned_1.fq
# REQUIRES: trimAdaptersSingle.py, trimPolyATSingle.py, FastQC (if fastqc command previously run)
#
# PAIRED-END READS:
# INPUT: $SAMPLE/raw/unaligned_1.fq and $SAMPLE/raw/unaligned_2.fq
# OUTPUT: $SAMPLE/trimAT/unaligned_1.fq and $SAMPLE/trimAT/unaligned_2.fq, $SAMPLE/trimAdapters.stats.txt, $SAMPLE/trimPolyAT.stats.txt
#         intermediate files $SAMPLE/trimAD/unaligned_1.fq and $SAMPLE/trimAD/unaligned_2.fq
# REQUIRES: trimReads.py, FastQC (if fastqc command previously run)
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_TRIM="Usage: `basename $0` trim OPTIONS sampleID    --  trim adapter and poly A/T contamination\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_TRIM="Usage: `basename $0` trim [-i inputDir] [-m minLen] [-se] sampleID\n"
ngsHelp_TRIM+="Input:\n\t$REPO_LOCATION/INPUTDIR/contaminants.fa (contaminants file)\n\tsampleID/orig/unaligned_1.fq\n\tsampleID/orig/unaligned_2.fq (paired-end reads)\n"
ngsHelp_TRIM+="Output:\n\tsampleID/INPUTDIR/unaligned_1.fq\n\tsampleID/INPUTDIR/unaligned_2.fq (paired-end reads)\n\tsampleID/INPUTDIR/stats.txt\n\tsampleID/INPUTDIR/contaminants.fa (contaminants file)\n"
ngsHelp_TRIM+="Requires:\n\ttrimReads.py\n\tFastQC (if fastqc command previously run)\n"
ngsHelp_TRIM+="Options:\n"
ngsHelp_TRIM+="\t-i inputDir - location of source files (default: orig).\n"
ngsHelp_TRIM+="\t-m minLen - Minimum size of trimmed read. If trimmed beyond minLen, then read is discarded. If read is paired then read is replaced with N's, unless both reads in pair are smaller than minLen in which case the pair is discarded. (default: 20).\n"
ngsHelp_TRIM+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_TRIM+="Runs trimReads.py to trim data. Trimmed data is placed in 'sampleID/INPUTDIR'."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# TRIM args: -se (optional), sampleID
##########################################################################################

ngsArgs_TRIM() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	fi

	# default value
	INPDIR="orig"
	MINLEN="20"

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) INPDIR=$2
				shift; shift;
				;;
			-p) MINLEN=$2
				shift; shift;
				;;
			-se) SE=true
				shift;
				;;
			-*) printf "Illegal option: '%s'\n" "$1"
				printHelp $COMMAND
				exit 0
				;;
 			*) break ;;
		esac
	done
	
	SAMPLE=$1
}

##########################################################################################
# RUNNING COMMAND ACTION
# TRIM command
##########################################################################################

ngsCmd_TRIM() {
	if $SE; then prnCmd "# BEGIN: TRIMMING SINGLE-END"
	else prnCmd "# BEGIN: TRIMMING PAIRED-END"; fi
		
	# make relevant directory
	if ! $DEBUG; then 
		prnCmd "mkdir $SAMPLE/trim"
		if [ ! -d $SAMPLE/trim ]; then mkdir $SAMPLE/trim; fi
	fi
	
	if $SE; then
		# single-end
		prnCmd "trimReads.py -m $MINLEN -rN -rAT 26 -c $REPO_LOCATION/trim/contaminants.fa -f $SAMPLE/$INPDIR/unaligned_1.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt"
		if ! $DEBUG; then 
			trimReads.py -m $MINLEN -rN -rAT 26 -c $REPO_LOCATION/trim/contaminants.fa -f $SAMPLE/$INPDIR/unaligned_1.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt
		fi
		
	else
		# paired-end
		prnCmd "trimReads.py -p -m $MINLEN -rN -rAT 26 -c $REPO_LOCATION/trim/contaminants.fa -f $SAMPLE/$INPDIR/unaligned_1.fq -r $SAMPLE/$INPDIR/unaligned_2.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt"
		if ! $DEBUG; then 
			trimReads.py -p -m $MINLEN -rN -rAT 26 -c $REPO_LOCATION/trim/contaminants.fa -f $SAMPLE/$INPDIR/unaligned_1.fq -r $SAMPLE/orig/unaligned_2.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt
		fi
	fi
	
	# if we ran fastqc on orig data (ie $SAMPLE/fastqc exists), then
	# also run on trimmed data. Put this fastqc output in separate
	# directory so it doesn't squash the output from orig
	if [ -d $SAMPLE/fastqc ]; then 
		if [ ! -d $SAMPLE/fastqc.trim ]; then 
			prnCmd "mkdir $SAMPLE/fastqc.trim"
			if ! $DEBUG; then mkdir $SAMPLE/fastqc.trim; fi
		fi
		prnCmd "fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trim/unaligned_1.fq"
		if ! $DEBUG; then 
			fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trim/unaligned_1.fq
			
			# do some cleanup of the output files
			prnCmd "mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/."
			mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/.
			
			prnCmd "rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc"
			rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc
		fi
	fi
	
	if $SE; then prnCmd "# FINISHED: TRIMMING SINGLE-END"
	else prnCmd "# FINISHED: TRIMMING PAIRED-END"; fi
}

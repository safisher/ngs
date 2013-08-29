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
# INPUT: $SAMPLE/orig/unaligned_1.fq
# OUTPUT: $SAMPLE/trim/unaligned_1.fq, $SAMPLE/trim/stats.txt
# REQUIRES: trimReads.py, FastQC (if fastqc command previously run)
#
# PAIRED-END READS:
# INPUT: $SAMPLE/orig/unaligned_1.fq and $SAMPLE/orig/unaligned_2.fq
# OUTPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq, $SAMPLE/trim/stats.txt
# REQUIRES: trimReads.py, FastQC (if fastqc command previously run)
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_TRIM="Usage: `basename $0` trim OPTIONS sampleID    --  trim adapter and poly A/T contamination\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_TRIM="Usage: `basename $0` trim [-i inputDir] [-c contaminantsFile] [-m minLen] [-se] sampleID\n"
ngsHelp_TRIM+="Input:\n\t$REPO_LOCATION/trim/contaminants.fa (file containing contaminants)\n\tsampleID/inputDir/unaligned_1.fq\n\tsampleID/inputDir/unaligned_2.fq (paired-end reads)\n"
ngsHelp_TRIM+="Output:\n\tsampleID/trim/unaligned_1.fq\n\tsampleID/trim/unaligned_2.fq (paired-end reads)\n\tsampleID/trim/stats.txt\n\tsampleID/trim/contaminants.fa (contaminants file)\n"
ngsHelp_TRIM+="Requires:\n\ttrimReads.py ( https://github.com/safisher/ngs )\n"
ngsHelp_TRIM+="Options:\n"
ngsHelp_TRIM+="\t-i inputDir - location of source files (default: orig).\n"
ngsHelp_TRIM+="\t-c contaminantsFile - file containing contaminants to be trimmed (default: $REPO_LOCATION/trim/contaminants.fa).\n"
ngsHelp_TRIM+="\t-m minLen - Minimum size of trimmed read. If trimmed beyond minLen, then read is discarded. If read is paired then read is replaced with N's, unless both reads in pair are smaller than minLen in which case the pair is discarded. (default: 20).\n"
ngsHelp_TRIM+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_TRIM+="Runs trimReads.py to trim data. Trimmed data is placed in 'sampleID/trim'. The contaminants file that was used is copied into the trim directory for future reference."

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_TRIM_INP_DIR="orig"
ngsLocal_TRIM_CONTAMINANTS_FILE="$REPO_LOCATION/trim/contaminants.fa"
ngsLocal_TRIM_MINLEN="20"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# TRIM args: -se (optional), sampleID
##########################################################################################

ngsArgs_TRIM() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_TRIM_INP_DIR=$2
				shift; shift;
				;;
			-c) ngsLocal_TRIM_CONTAMINANTS_FILE=$2
				shift; shift;
				;;
			-p) ngsLocal_TRIM_MINLEN=$2
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
	prnCmd "mkdir $SAMPLE/trim"
	if ! $DEBUG; then 
		if [ ! -d $SAMPLE/trim ]; then mkdir $SAMPLE/trim; fi
	fi
	
	# print version info in journal file
	prnCmd "# trimReads.py version"
	if ! $DEBUG; then prnCmd "# `trimReads.py -v`"; fi

	if $SE; then
		# single-end
		prnCmd "trimReads.py -m $ngsLocal_TRIM_MINLEN -rN -rAT 26 -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt"
		if ! $DEBUG; then 
			trimReads.py -m $ngsLocal_TRIM_MINLEN -rN -rAT 26 -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt
		fi
		
	else
		# paired-end
		prnCmd "trimReads.py -p -m $ngsLocal_TRIM_MINLEN -rN -rAT 26 -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -r $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_2.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt"
		if ! $DEBUG; then 
			trimReads.py -p -m $ngsLocal_TRIM_MINLEN -rN -rAT 26 -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -r $SAMPLE/orig/unaligned_2.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/stats.txt
		fi
	fi
	
	# copy contaminants files into trim directory for future reference
	prnCmd "cp $ngsLocal_TRIM_CONTAMINANTS_FILE $SAMPLE/trim/."
	if ! $DEBUG; then
		cp $ngsLocal_TRIM_CONTAMINANTS_FILE $SAMPLE/trim/.
	fi

	if $SE; then prnCmd "# FINISHED: TRIMMING SINGLE-END"
	else prnCmd "# FINISHED: TRIMMING PAIRED-END"; fi
}

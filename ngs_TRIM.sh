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
# INPUT: $SAMPLE/init/unaligned_1.fq
# OUTPUT: $SAMPLE/trim/unaligned_1.fq, $SAMPLE/trim/$SAMPLE.trim.stats.txt
# REQUIRES: trimReads.py, FastQC (if fastqc command previously run)
#
# PAIRED-END READS:
# INPUT: $SAMPLE/init/unaligned_1.fq and $SAMPLE/init/unaligned_2.fq
# OUTPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq, $SAMPLE/trim/$SAMPLE.trim.stats.txt
# REQUIRES: trimReads.py, FastQC (if fastqc command previously run)
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` trim OPTIONS sampleID    --  trim adapter and poly A/T contamination\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_TRIM() {
	echo -e "Usage: `basename $0` trim [-i inputDir] [-c contaminantsFile] [-m minLen] [-kN] [-rAT numBases] [-se] sampleID"
	echo -e "Input:\n\t$REPO_LOCATION/trim/contaminants.fa (file containing contaminants)\n\tsampleID/inputDir/unaligned_1.fq\n\tsampleID/inputDir/unaligned_2.fq (paired-end reads)"
	echo -e "Output:\n\tsampleID/trim/unaligned_1.fq\n\tsampleID/trim/unaligned_2.fq (paired-end reads)\n\tsampleID/trim/sampleID.trim.stats.txt\n\tsampleID/trim/contaminants.fa (contaminants file)"
	echo -e "Requires:\n\ttrimReads.py ( https://github.com/safisher/ngs )"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source files (default: init)."
	echo -e "\t-c contaminantsFile - file containing contaminants to be trimmed (default: $REPO_LOCATION/trim/contaminants.fa)."
	echo -e "\t-m minLen - Minimum size of trimmed read. If trimmed beyond minLen, then read is discarded. If read is paired then read is replaced with N's, unless both reads in pair are smaller than minLen in which case the pair is discarded. (default: 20)."
	echo -e "\t-kN - do not trim N's from either end of the reads (default: remove N's)."
	echo -e "\t-rAT numBases - number of polyA/T bases to trim (default: 26). To disable polyA/T trimming, use '0' (eg '-rAT 0')."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "Runs trimReads.py to trim data. Trimmed data is placed in 'sampleID/trim'. The contaminants file that was used is copied into the trim directory for future reference."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_TRIM_INP_DIR="init"
ngsLocal_TRIM_CONTAMINANTS_FILE="$REPO_LOCATION/trim/contaminants.fa"
ngsLocal_TRIM_MINLEN="20"
ngsLocal_TRIM_POLYAT_BASES="26"
ngsLocal_TRIM_KN=false

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# TRIM args: -se (optional), sampleID
##########################################################################################

ngsArgs_TRIM() {
	if [ $# -lt 1 ]; then printHelp "TRIM"; fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_TRIM_INP_DIR=$2
				shift; shift;
				;;
			-c) ngsLocal_TRIM_CONTAMINANTS_FILE=$2
				shift; shift;
				;;
			-m) ngsLocal_TRIM_MINLEN=$2
				shift; shift;
				;;
			-kN) ngsLocal_TRIM_KN=true
				shift;
				;;
			-rAT) ngsLocal_TRIM_POLYAT_BASES=$2
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
	if [ ! -d $SAMPLE/trim ]; then 
		prnCmd "mkdir $SAMPLE/trim"
		if ! $DEBUG; then mkdir $SAMPLE/trim; fi
	fi
	
    # print version info in $SAMPLE directory
	prnCmd "# trimReads.py version: trimReads.py -v 2>&1"
	if ! $DEBUG; then 
		ver=$(trimReads.py -v 2>&1)
		prnVersion "trim" "program\tversion" "trimReads.py\t$ver"
	fi

	# set argument for removing of N's from both ends of the reads
	ngsLocal_TRIM_RN="-rN"
	if $ngsLocal_TRIM_KN; then
		ngsLocal_TRIM_RN=""
	fi

	# set polyA/T trimming argument
	ngsLocal_TRIM_POLYAT=""
	if [[ $ngsLocal_TRIM_POLYAT_BASES -gt 0 ]]; then
		ngsLocal_TRIM_POLYAT="-rAT $ngsLocal_TRIM_POLYAT_BASES"
	fi

	if $SE; then
		# single-end
		prnCmd "trimReads.py -m $ngsLocal_TRIM_MINLEN $ngsLocal_TRIM_RN $ngsLocal_TRIM_POLYAT -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/$SAMPLE.trim.stats.txt"
		if ! $DEBUG; then 
			trimReads.py -m $ngsLocal_TRIM_MINLEN $ngsLocal_TRIM_RN $ngsLocal_TRIM_POLYAT -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/$SAMPLE.trim.stats.txt
		fi
		
	else
		# paired-end
		prnCmd "trimReads.py -p -m $ngsLocal_TRIM_MINLEN $ngsLocal_TRIM_RN $ngsLocal_TRIM_POLYAT -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -r $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_2.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/$SAMPLE.trim.stats.txt"
		if ! $DEBUG; then 
			trimReads.py -p -m $ngsLocal_TRIM_MINLEN $ngsLocal_TRIM_RN $ngsLocal_TRIM_POLYAT -c $ngsLocal_TRIM_CONTAMINANTS_FILE -f $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq -r $SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_2.fq -o $SAMPLE/trim/unaligned > $SAMPLE/trim/$SAMPLE.trim.stats.txt
		fi
	fi
	
	# copy contaminants files into trim directory for future reference
	prnCmd "cp $ngsLocal_TRIM_CONTAMINANTS_FILE $SAMPLE/trim/."
	if ! $DEBUG; then
		cp $ngsLocal_TRIM_CONTAMINANTS_FILE $SAMPLE/trim/.
	fi

	# run error checking
	if ! $DEBUG; then ngsErrorChk_TRIM $@; fi

	if $SE; then prnCmd "# FINISHED: TRIMMING SINGLE-END"
	else prnCmd "# FINISHED: TRIMMING PAIRED-END"; fi
}

##########################################################################################
# ERROR CHECKING. Confirm that the input files are not empty and if PE
# then make sure the files have the same number of lines.
##########################################################################################

ngsErrorChk_TRIM() {
	prnCmd "# TRIM ERROR CHECKING: RUNNING"

	inputFile_1="$SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_1.fq"
	outputFile_1="$SAMPLE/trim/unaligned_1.fq"

	if $SE; then
		# make sure expected output file exists and is not empty
		if [ ! -s $outputFile_1 ]; then
			errorMsg="Expected TRIM output file does not exist.\n"
			errorMsg+="\tinput file: $inputFile_1\n"
			errorMsg+="\toutput file: $outputFile_1\n"
			prnError "$errorMsg"
		fi

	else
		# paired-end
		inputFile_2="$SAMPLE/$ngsLocal_TRIM_INP_DIR/unaligned_2.fq"
		outputFile_2="$SAMPLE/trim/unaligned_2.fq"

		# make sure expected output files exists
		if [[ ! -s $outputFile_1 || ! -s $outputFile_2 ]]; then
			errorMsg="Error with TRIM output files (don't exist or are empty).\n"
			errorMsg+="\tinput file: $inputFile_1\n"
			errorMsg+="\toutput file: $outputFile_1\n\n"
			errorMsg+="\tinput file: $inputFile_2\n"
			errorMsg+="\toutput file: $outputFile_2\n"
			prnError "$errorMsg"
		fi

		# compute number of lines in output files
		mate1=`wc -l $outputFile_1 | awk '{print $1}'`
		mate2=`wc -l $outputFile_2 | awk '{print $1}'`

		# make sure read files have same number of lines.
		if [ "$mate1" -ne "$mate2" ]; then
			errorMsg="Trimmed output files do not have the same number of lines.\n"
			errorMsg+="\tnum lines in first read file: $mate1\n"
			errorMsg+="\tinput file: $inputFile_1\n"
			errorMsg+="\toutput file: $outputFile_1\n\n"
			errorMsg+="\tnum lines in second read file: $mate2\n"
			errorMsg+="\tinput file: $inputFile_2\n"
			errorMsg+="\toutput file: $outputFile_2\n"
			prnError "$errorMsg"
		fi
	fi

	prnCmd "# TRIM ERROR CHECKING: DONE"
}

##########################################################################################
# PRINT STATS. Prints a tab-delimited list stats of interest.
##########################################################################################

ngsStats_TRIM() {
	if [ $# -ne 1 ]; then
		prnError "Incorrect number of parameters for ngsStats_TRIM()."
	fi

	case $1 in
		header)
			# the second to the last line of the stats file is a tab-delimited lists of headers
			echo `tail -2 $SAMPLE/trim/$SAMPLE.trim.stats.txt | head -1`
			;;

		values)
			# the last line of the stats file is a tab-delimited lists of values
			echo `tail -1 $SAMPLE/trim/$SAMPLE.trim.stats.txt`
			;;

		keyvalue)
			# output key:value pair of stats

			# the bash IFS variable dictates the word delimiting which is " \t\n" 
			# by default. We want to only delimite by tabs for the case here.
			local IFS=$'\t'
			
			# the last two lines of the stats.txt file are tab-delimited lists of headers and values
			declare -a header=($(tail -2 $SAMPLE/trim/$SAMPLE.trim.stats.txt | head -1))
			declare -a values=($(tail -1 $SAMPLE/trim/$SAMPLE.trim.stats.txt))
			
			# output a tab-delimited, key:value list
			numFields=${#header[@]}
			for ((i=0; i<$numFields-1; ++i)); do
				echo -en "${header[$i]}:${values[$i]}\t"
			done
			echo "${header[$numFields-1]}:${values[$numFields-1]}"
			;;

		*) 
			# incorrect argument
			prnError "Invalid parameter for ngsStats_TRIM() (got $1, expected: 'header|values')."
			;;
	esac
}

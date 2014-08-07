#!/bin/bash

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

##########################################################################################
# SINGLE-END READS
# INPUT: $SAMPLE/init/unaligned_1.fq
# OUTPUT: $SAMPLE/rmdup/*
#
# PAIRED-END READS
# INPUT: $SAMPLE/init/unaligned_1.fq and $SAMPLE/init/unaligned_2.fq
# OUTPUT: $SAMPLE/rmdup/*
#
# REQUIRES: removeDuplicates.py
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` rmdup OPTIONS sampleID    --  remove duplicate reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_RMDUP() {
	echo -e "Usage:\n\t`basename $0` rmdup [-i inputDir] [-se] sampleID"
	echo -e "Input:\n\tsampleID/inputDir/unaligned_1.fq\n\tsampleID/inputDir/unaligned_2.fq (paired-end reads)"
	echo -e "Output:\n\tsampleID/rmdup/unaligned_1.fq\n\tsampleID/rmdup/unaligned_1.fq\n \tsampleID/rmdup/sampleID.rmdup.stats.txt"
	echo -e "Requires:\n\tremoveDuplicates.py"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source files (default: init)."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "Remove duplicate reads. Reads are considered duplicates if they exactly match. For paired-end reads, the mate pairs both must exactly match to be considered duplicates. This is very RAM intensive, requiring RAM amounts up to three times the input file size (e.g. if your fastq files total 20GB then up to 60GB RAM may be used when removing duplicates)."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_RMDUP_INP_DIR="init"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# RMDUP args: -se (optional), sampleID
##########################################################################################

ngsArgs_RMDUP() {
	if [ $# -lt 1 ]; then printHelp "RMDUP"; fi
	
	SAMPLE=$1
	
	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_RMDUP_INP_DIR=$2
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
# Remove duplicate reads.
##########################################################################################

ngsCmd_RMDUP() {
	if $SE; then prnCmd "# BEGIN: REMOVE DUPLICATES SINGLE-END"
	else prnCmd "# BEGIN: REMOVE DUPLICATES PAIRED-END"; fi
	
	# make relevant directory
	if [ ! -d $SAMPLE/rmdup ]; then 
		prnCmd "mkdir $SAMPLE/rmdup"
		if ! $DEBUG; then mkdir $SAMPLE/rmdup; fi
	fi
	
	prnCmd "# removeDuplicates.py version: removeDuplicates.py -v"
	if ! $DEBUG; then 
		ver=$(removeDuplicates.py -v 2>&1)
		prnVersion "rmdup" "program\tversion" "removeDuplicates.py\t$ver"
	fi

	# Convert init/fastq files into single fasta file (raw.fa)
	if $SE; then
		prnCmd "removeDuplicates.py -f $SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_1.fq -o $SAMPLE/rmdup/unaligned > $SAMPLE/rmdup/$SAMPLE.rmdup.stats.txt"
		if ! $DEBUG; then 
			removeDuplicates.py -f $SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_1.fq -o $SAMPLE/rmdup/unaligned > $SAMPLE/rmdup/$SAMPLE.rmdup.stats.txt
		fi
	else
		prnCmd "removeDuplicates.py -f $SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_1.fq -r $SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_2.fq -o $SAMPLE/rmdup/unaligned > $SAMPLE/rmdup/$SAMPLE.rmdup.stats.txt"
		if ! $DEBUG; then 
			removeDuplicates.py -f $SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_1.fq -r $SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_2.fq -o $SAMPLE/rmdup/unaligned > $SAMPLE/rmdup/$SAMPLE.rmdup.stats.txt
		fi
	fi
	
	if $SE; then prnCmd "# FINISHED: REMOVE DUPLICATES SINGLE-END"
	else prnCmd "# FINISHED: REMOVE DUPLICATES PAIRED-END"; fi
}

##########################################################################################
# ERROR CHECKING. Make sure output file exists and contains species counts.
##########################################################################################

ngsErrorChk_RMDUP() {
	prnCmd "# RMDUP ERROR CHECKING: RUNNING"

	inputFile_1="$SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_1.fq"
	inputFile_2="$SAMPLE/$ngsLocal_RMDUP_INP_DIR/unaligned_2.fq"
	outputFile_1="$SAMPLE/rmdup/unaligned_1.fq"
	outputFile_2="$SAMPLE/rmdup/unaligned_2.fq"
	outputFile_3="$SAMPLE/rmdup/$SAMPLE.rmdup.stats.txt"

	# make sure expected output files exists
	if [[ ! -s $outputFile_1 || ! -s $outputFile_2 || ! -s $outputFile_3 ]]; then
		errorMsg="Error with RMDUP output files (don't exist or are empty).\n"
		errorMsg+="\tinput file: $inputFile_1\n"
		if ! $SE; then errorMsg+="\tinput file: $inputFile_2\n"; fi
		errorMsg+="\toutput file: $outputFile_1\n"
		if ! $SE; then errorMsg+="\toutput file: $outputFile_2\n"; fi
		errorMsg+="\tstats file: $outputFile_3\n"
		prnError "$errorMsg"
	fi

	# compute number of lines in stats file
	counts=`wc -l $outputFile_3 | awk '{print $1}'`

	# if counts file has less than 3 lines, then RMDUP didn't work
	if [ "$counts" -lt "3" ]; then
		errorMsg="RMDUP failed to run properly and there are no stats.\n"
		errorMsg+="\tinput file: $inputFile_1\n"
		if ! $SE; then errorMsg+="\tinput file: $inputFile_2\n"; fi
		errorMsg+="\toutput file: $outputFile_1\n"
		if ! $SE; then errorMsg+="\toutput file: $outputFile_2\n"; fi
		errorMsg+="\tstats file: $outputFile_3\n"
		prnError "$errorMsg"
	fi

	prnCmd "# RMDUP ERROR CHECKING: DONE"
}

##########################################################################################
# PRINT STATS. Prints a tab-delimited list stats of interest.
##########################################################################################

ngsStats_RMDUP() {
	if [ $# -ne 1 ]; then
		prnError "Incorrect number of parameters for ngsStats_RMDUP()."
	fi
	
	local statsFile="$SAMPLE/rmdup/$SAMPLE.rmdup.stats.txt"

	# the second to the last line of the stats file is a tab-delimited lists of headers
	local header=$(tail -2 $statsFile | head -1)
	# the last line of the stats file is a tab-delimited lists of values
	local values=$(tail -1 $statsFile)

	case $1 in
		header)
			echo $header
			;;

		values) 
			echo $values
			;;

		keyvalue) 
			# output key:value pair of stats

			# the bash IFS variable dictates the word delimiting which is " \t\n" 
			# by default. We want to only delimite by tabs for the case here.
			local IFS=$'\t'

			# convert tab-delimited header/values variables to array
			declare -a headerArray=($header)
			declare -a valuesArray=($values)
			
			# output a tab-delimited, key:value list
			numFields=${#headerArray[@]}
			for ((i=0; i<$numFields-1; ++i)); do
				echo -en "${headerArray[$i]}:${valuesArray[$i]}\t"
			done
			echo "${headerArray[$numFields-1]}:${valuesArray[$numFields-1]}"
			;;

		*) 
			# incorrect argument
			prnError "Invalid parameter for ngsStats_RMDUP() (got $1, expected: 'header|values')."
			;;
	esac
}

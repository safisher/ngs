#!/bin/bash

# Copyright (c) 2012,2013,2016 Stephen Fisher, Jamie Shallcross, Junhyong Kim 
# University of Pennsylvania.  All Rights Reserved.
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
# INPUT: list input files here
# OUTPUT: list output files here
# REQUIRES: list any external programs required to complete COMMAND function
##########################################################################################

##########################################################################################
# USAGE
# ngsUsage_BARCODE should be a single line that ends in a "\n"
##########################################################################################

NGS_USAGE+="Return all reads that contain the specified feature, allowing for one mismatch. The k-mer is removed from the read if the -c flag is present.\n"

##########################################################################################
# HELP TEXT
# ngsHelp_BARCODE should contain expanded help
##########################################################################################

ngsHelp_BARCODE() {
	echo -e "Usage: `basename $0` barcode [-i inputDir] [-t numProc] [-p] [-c contaminantsFile] [-m minLen] [-q phredThreshold] [-rN] [-rAT numBases] [-se] -b barcode sampleID"
	echo -e "Input:\n\t$REPO_LOCATION/trim/contaminants.fa (file containing contaminants)\n\tsampleID/inputDir/unaligned_1.fq\n\tsampleID/inputDir/unaligned_2.fq (paired-end reads)"
	echo -e "Output:\n\tsampleID/trim/unaligned_1.fq\n\tsampleID/trim/unaligned_2.fq (paired-end reads)\n\tsampleID/trim/sampleID.trim.stats.txt\n\tsampleID/trim/contaminants.fa (contaminants file)"
	echo -e "Requires:\n\ttrimReads.py ( https://github.com/safisher/ngs )"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source files (default: init)."
	echo -e "\t-t numProc - maximum number of cpu to use."
	echo -e "\t-c - Cut reads at the barcode point. Read 1 will be cut from the barcode to the 3' end, read 2 will be cut from the barcode to the 5' end."
	echo -e "\t-m - Mask the barcode portion of reads. Sequences on both sides of the barcode will remain unchanged."
	echo -e "\t-l length - length of reads"
	echo -e "\t-pre prefix - sequence to be inserted before barcode sequence (default: none)"
	echo -e "\t-suf suffix - sequence to be inserted after barcode sequence (default: none)"
	echo -e "\t-m [AND|OR] - if a both a prefix and a suffix is given, whether to look for reads with [prefix][barcode][suffix] (AND), or reads with either [prefix][barcode] or [barcode][suffix] (OR)."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "\t-b barcode - barcode to be trimmed."
	echo -e "Runs kmerFind.py to select reads with a specific barcode sequence, optionally trimming with -c. Selected/trimmed data is placed in 'sampleID/barcode.trim'."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

# put local variables here using module naming convention: ngsLocal_BARCODE_VARIABLENAME

ngsLocal_BARCODE_INP_DIR=$SAMPLE/init
ngsLocal_BARCODE_CUT_READS=""
#ngsLocal_BARCODE_R1_THRESH=45
#ngsLocal_BARCODE_R2_THRESH=10
ngsLocal_BARCODE_R1_THRESH=0 # no threshold, finds anywhere in read
ngsLocal_BARCODE_R2_THRESH=0 # no threshold, finds anywhere in read
ngsLocal_BARCODE_PREFIX=""
ngsLocal_BARCODE_SUFFIX=""
ngsLocal_BARCODE_MODE=""

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# BARCODE args
##########################################################################################

ngsArgs_BARCODE() {
	# process arguments here
	if [ $# -lt 1 ]; then printHelp "BARCODE"; fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_BARCODE_INP_DIR=$2
				shift; shift;
				;;
			-b) ngsLocal_BARCODE_BARCODE=$2
				shift; shift;
				;;
			-c) ngsLocal_BARCODE_CUT_READS="-c"
				shift;
				;;
			-m) ngsLocal_BARCODE_CUT_READS="-m"
				shift;
				;;
			-l) READ_LENGTH=$2
				shift; shift;
				;;
			-se) SE=true
				shift;
				;;
			-pre) ngsLocal_BARCODE_PREFIX="-prefix $2"
				shift; shift;
				;;
			-suf) ngsLocal_BARCODE_SUFFIX="-suffix $2"
				shift; shift;
				;;
			-m) ngsLocal_BARCODE_MODE="-psMode $2"
				shift; shift;
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
# BARCODE command
##########################################################################################

ngsCmd_BARCODE() {
	# 1. create subdirectory, if necessary
	# 2. output version information using prnVersion() (see ngs.sh)
	# 3. do stuff here
	# 4. call ngsErrorChk_BARCODE

	if $SE; then prnCmd "# begin: kmer finding single-end"
	else prnCmd "# begin: kmer finding paired-end"; fi
    
	# make relevant directory
	if [ ! -d $SAMPLE/barcode.trim ]; then 
		prnCmd "mkdir $SAMPLE/barcode.trim"
		if ! $DEBUG; then mkdir $SAMPLE/barcode.trim; fi
	fi
    
	# uncomment if using thresholds
	#	if [[ $READ_LENGTH -ge 100 ]]; then
	#	    ngsLocal_BARCODE_R1_THRESH=75
	#	fi
    	
	# print version info in $sample directory
	prnCmd "# kmerFind.py version: kmerFind.py -v 2>&1"
	if ! $DEBUG; then 
		local prnSE=0
		if $SE; then prnSE=1; fi
		local prnCut=0;
		if [[ $ngsLocal_BARCODE_CUT_READS == "-c" ]]; then prnCut=1; fi
		ver=$(kmerFind.py -v 2>&1)
		prnVersion "barcode.trim" \
		"program\tversion\tprefix\tbarcode\tsuffix\tcut\tSE\tR1_threshold\tR2_threshold" \
		"kmerFind.py\t$ver\t$ngsLocal_BARCODE_PREFIX\t$ngsLocal_BARCODE_BARCODE\t$ngsLocal_BARCODE_SUFFIX\t$prnCut\t$prnSE\t$ngsLocal_BARCODE_R1_THRESH\t$ngsLocal_BARCODE_R2_THRESH"
	fi
    

	if $SE; then
	    prnCmd "kmerFind.py -TF $ngsLocal_BARCODE_R1_THRESH -b $ngsLocal_BARCODE_BARCODE $ngsLocal_BARCODE_PREFIX $ngsLocal_BARCODE_SUFFIX $ngsLocal_BARCODE_MODE -t 6 $ngsLocal_BARCODE_CUT_READS -f $SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_1.fq -o $SAMPLE/barcode.trim/unaligned > $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt" 
	    if ! $DEBUG; then 
		kmerFind.py -TF $ngsLocal_BARCODE_R1_THRESH -b $ngsLocal_BARCODE_BARCODE $ngsLocal_BARCODE_PREFIX $ngsLocal_BARCODE_SUFFIX $ngsLocal_BARCODE_MODE -t 6 $ngsLocal_BARCODE_CUT_READS -f $SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_1.fq -o $SAMPLE/barcode.trim/unaligned > $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt
	    fi
	else
	    prnCmd "kmerFind.py -TF $ngsLocal_BARCODE_R1_THRESH -TR $ngsLocal_BARCODE_R2_THRESH -b $ngsLocal_BARCODE_BARCODE $ngsLocal_BARCODE_PREFIX $ngsLocal_BARCODE_SUFFIX $ngsLocal_BARCODE_MODE -t 6 $ngsLocal_BARCODE_CUT_READS -f $SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_1.fq -r $SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_2.fq -o $SAMPLE/barcode.trim/ > $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt" 
	    if ! $DEBUG; then 
		kmerFind.py -TF $ngsLocal_BARCODE_R1_THRESH -TR $ngsLocal_BARCODE_R2_THRESH -b $ngsLocal_BARCODE_BARCODE $ngsLocal_BARCODE_PREFIX $ngsLocal_BARCODE_SUFFIX $ngsLocal_BARCODE_MODE -t 6 $ngsLocal_BARCODE_CUT_READS -f $SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_1.fq -r $SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_2.fq -o $SAMPLE/barcode.trim/ > $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt
	    fi
	fi
	if [[ -f $SAMPLE/barcode.trim/unaligned_1.hist.pdf ]]; then #we test this b/c if there are no mapped reads no histogram is written
	    prnCmd "mv $SAMPLE/barcode.trim/unaligned_1.hist.pdf $SAMPLE/barcode.trim/${SAMPLE}.hist_1.pdf"
	    if ! $DEBUG; then
		mv $SAMPLE/barcode.trim/unaligned_1.hist.pdf $SAMPLE/barcode.trim/${SAMPLE}.hist_1.pdf
	    fi
	fi
	if [[ -f $SAMPLE/barcode.trim/unaligned_2.hist.pdf ]]; then
	    prnCmd "mv $SAMPLE/barcode.trim/unaligned_2.hist.pdf $SAMPLE/barcode.trim/${SAMPLE}.hist_2.pdf"
	    if ! $DEBUG; then
		mv $SAMPLE/barcode.trim/unaligned_2.hist.pdf $SAMPLE/barcode.trim/${SAMPLE}.hist_2.pdf
	    fi
	fi

	prnCmd "# FINISHED: BARCODE"
}

##########################################################################################
# ERROR CHECKING
##########################################################################################

ngsErrorChk_BARCODE() {
    # check to make sure output files are valid (eg they exist and are not empty)
    prnCmd "# BARCODE ERROR CHECKING: RUNNING"

    inputFile_1="$SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_1.fq"
    outputFile_1="$SAMPLE/barcode.trim/unaligned_1.fq"

    if $SE; then
	# make sure expected output file exists and is not empty
	if [ ! -s $outputFile_1 ]; then
	    errorMsg="Expected BARCODE output file does not exist.\n"
	    errorMsg+="\tinput file: $inputFile_1\n"
	    errorMsg+="\toutput file: $outputFile_1\n"
	    prnError "$errorMsg"
	fi

    else
	# paired-end
	inputFile_2="$SAMPLE/$ngsLocal_BARCODE_INP_DIR/unaligned_2.fq"
	outputFile_2="$SAMPLE/barcode.trim/unaligned_2.fq"

	# make sure expected output files exists
	if [[ ! -s $outputFile_1 || ! -s $outputFile_2 ]]; then
	    errorMsg="Error with BARCODE output files (don't exist or are empty).\n"
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
	    errorMsg="barcodemed output files do not have the same number of lines.\n"
	    errorMsg+="\tnum lines in first read file: $mate1\n"
	    errorMsg+="\tinput file: $inputFile_1\n"
	    errorMsg+="\toutput file: $outputFile_1\n\n"
	    errorMsg+="\tnum lines in second read file: $mate2\n"
	    errorMsg+="\tinput file: $inputFile_2\n"
	    errorMsg+="\toutput file: $outputFile_2\n"
	    prnError "$errorMsg"
	fi
    fi

    prnCmd "# BARCODE ERROR CHECKING: DONE"
}

##########################################################################################
# PRINT STATS
##########################################################################################

ngsStats_BARCODE() {
    # expect one argument that is either "header" or "values" and output tab-delimited list stats of interest.


    if [ $# -ne 1 ]; then
	prnError "Incorrect number of parameters for ngsStats_BARCODE()."
    fi

    case $1 in
	header)
	    # the second to the last line of the stats file is a tab-delimited lists of headers
	    echo `tail -2 $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt | head -1`
	    ;;

	values)
	    # the last line of the stats file is a tab-delimited lists of values
	    echo `tail -1 $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt`
	    ;;

	keyvalue)
	    # output key:value pair of stats

	    # the bash IFS variable dictates the word delimiting which is " \t\n" 
	    # by default. We want to only delimite by tabs for the case here.
	    local IFS=$'\t'

	    # the last two lines of the stats.txt file are tab-delimited lists of headers and values
	    declare -a header=($(tail -2 $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt | head -1))
	    declare -a values=($(tail -1 $SAMPLE/barcode.trim/$SAMPLE.barcode.stats.txt))

	    # output a tab-delimited, key:value list
	    numFields=${#header[@]}
	    for ((i=0; i<$numFields-1; ++i)); do
		echo -en "${header[$i]}:${values[$i]}\t"
	    done
	    echo "${header[$numFields-1]}:${values[$numFields-1]}"
	    ;;

	*) 
	    # incorrect argument
	    prnError "Invalid parameter for ngsStats_BARCODE() (got $1, expected: 'header|values')."
	    ;;
    esac

    wait
}

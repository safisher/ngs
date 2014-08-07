#!/bin/bash

# Copyright (c) 2012-2014, Stephen Fisher, Hoa Giang and Junhyong Kim, University of
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
# INPUT: list input files here
# OUTPUT: list output files here
# REQUIRES: list any external programs required to complete COMMAND function
##########################################################################################

# "SPAdes" is a place holder and should be replaced by the name of the command

##########################################################################################
# USAGE 
##########################################################################################

NGS_USAGE+="Usage: `basename $0` spades OPTIONS sampleID -- run SPADES on trimmed reads.\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_SPADES() {
	echo -e "Usage:\n\t`basename $0` spades [-i inputDir] [-m maxRAM] [-k kmers] -p numProc [-se] sampleID"
	echo -e "Input:\n\tsampleID/inputDir/notMapped_1.fq\n\tsampleID/inputDir/notMapped_2.fq (paired-end reads)"
	echo -e "Output:\n\tsampleID/spades/SampleID.fasta "
	echo -e "Requires:\n\tSPAdes ( http://bioinf.spbau.ru/spades )"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source files (default: bowtie)."
	echo -e "\t-p numProc - number of cpu to use."
	echo -e "\t-k kmers (default: 33,49,83) --> SPAdes"
	echo -e "\t-m maximum RAM (default: 250)"
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "Runs SPAdes using the reads that BOWTIE was not able to map (sampleID/bowtie). Output is stored in sampleID/spades."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_SPADES_INP_DIR="bowtie"
ngsLocal_SPADES_MAX_RAM=100
ngsLocal_SPADES_KMERS="33,49,83"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# SPADES args: -i value, -p value, -m value, -k value, -se (optional), sampleID
##########################################################################################

ngsArgs_SPADES() {
	if [ $# -lt 3 ]; then printHelp "SPADES"; fi
	
	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_SPADES_INP_DIR=$2
				shift; shift;
				;;
			-p) NUMCPU=$2
				shift; shift;
				;;
			-m) ngsLocal_SPADES_MAX_RAM=$2
				shift; shift;
				;;
			-k) ngsLocal_SPADES_KMERS=$2
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
# SPADES command, version 2.2.1
##########################################################################################

ngsCmd_SPADES() {
	if $SE; then prnCmd "# BEGIN: SPADES SINGLE-END ASSEMBLY"
	else prnCmd "# BEGIN: SPADES PAIRED-END ASSEMBLY"; fi
	
	# make relevant directory
	if [ ! -d $SAMPLE/spades ]; then 
		prnCmd "mkdir $SAMPLE/spades"
		if ! $DEBUG; then mkdir $SAMPLE/spades; fi
	fi
	
	# print version info in journal file
	prnCmd "# SPAdes version: (check SPAdes log file)"
	if ! $DEBUG; then 
		prnVersion "spades" "program\tversion" "spades.py\t2.2.1"
	fi

	if $SE; then
		# single-end
		prnCmd "spades.py -1 $SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_1.fq -t $NUMCPU -m $ngsLocal_SPADES_MAX_RAM -k $ngsLocal_SPADES_KMERS -n $SAMPLE -o $SAMPLE/spades"
		if ! $DEBUG; then 
			spades.py -1 $SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_1.fq -t $NUMCPU -m $ngsLocal_SPADES_MAX_RAM -k $ngsLocal_SPADES_KMERS -n $SAMPLE -o $SAMPLE/spades
		fi
	else
		# paired-end
		prnCmd "spades.py -1 $SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_1.fq -2 $SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_2.fq -t $NUMCPU -m $ngsLocal_SPADES_MAX_RAM -k $ngsLocal_SPADES_KMERS -n $SAMPLE -o $SAMPLE/spades"
		if ! $DEBUG; then 
			spades.py -1 $SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_1.fq -2 $SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_2.fq -t $NUMCPU -m $ngsLocal_SPADES_MAX_RAM -k $ngsLocal_SPADES_KMERS -n $SAMPLE -o $SAMPLE/spades
		fi
	fi

	# move the SPAdes output to be more accessible
	prnCmd "mv $SAMPLE/spades/$SAMPLE/contigs/$SAMPLE.fasta $SAMPLE/spades/$SAMPLE.fasta"
	if ! $DEBUG; then mv $SAMPLE/spades/$SAMPLE/contigs/$SAMPLE.fasta $SAMPLE/spades/$SAMPLE.fasta; fi

	# run error checking
	if ! $DEBUG; then ngsErrorChk_SPADES $@; fi
	
	if $SE; then prnCmd "# FINISHED: SPADES SINGLE-END ASSEMBLY"
	else prnCmd "# FINISHED: SPADES PAIRED-END ASSEMBLY"; fi
}

##########################################################################################
# ERROR CHECKING. Make sure output files exist and are not empty.
##########################################################################################

ngsErrorChk_SPADES() {
	prnCmd "# SPADES ERROR CHECKING: RUNNING"

	inputFile_1="$SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_1.fq"
	inputFile_2="$SAMPLE/$ngsLocal_SPADES_INP_DIR/notMapped_2.fq"
	outputFile="$SAMPLE/spades/$SAMPLE.fasta"

	# make sure expected output file exists and is not empty
	if [ ! -s $outputFile ]; then
		errorMsg="Expected output file does not exist or is empty.\n"
		errorMsg+="\tinput file: $inputFile_1\n"
		if ! $SE; then errorMsg+="\tinput file: $inputFile_2\n"; fi
		errorMsg+="\toutput file: $outputFile\n"
		prnError "$errorMsg"
	fi

	prnCmd "# SPADES ERROR CHECKING: DONE"
}

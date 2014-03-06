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

ngsUsage_SPAdes="Usage:\n\t`basename $0` SPAdes OPTIONS sampleID -- run SPAdes on trimmed reads.\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_SPAdes="Usage:\n\t`basename $0` SPAdes [-i inputDir] -p numProc -m maxRAM -k kmers [-se] sampleID\n"
ngsHelp_SPAdes+="Input:\n\tsampleID/inputDir/unaligned_1.fq\n\tsampleID/inputDir/unaligned_2.fq (paired-end reads)\n"
ngsHelp_SPAdes+="Output:\n\tsampleID/SPAdes/SampleID.fasta \n"
ngsHelp_SPAdes+="Requires:\n\tSPAdes ( http://bioinf.spbau.ru/spades )\n"
ngsHelp_SPAdes+="Options:\n"
ngsHelp_SPAdes+="\t-i inputDir - location of source files (default: trim).\n"
ngsHelp_SPAdes+="\t-p numProc - number of cpu to use.\n"
ngsHelp_SPAdes+="\t-k kmers (default: 33,49,83) --> SPAdes\n"
ngsHelp_SPAdes+="\t-m maximum RAM (default: 250)\n"
ngsHelp_SPAdes+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_SPAdes+="Runs SPAdes using the trimmed files from sampleID/trim. Output is stored in sampleID/SPAdes."

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_SPAdes_INP_DIR="trim"
ngsLocal_SPAdes_MAX_RAM=100
ngsLocal_SPAdes_KMERS="33,49,83"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# SPAdes args: -i value, -p value, -m value, -k value, -se (optional), sampleID
##########################################################################################

ngsArgs_SPAdes() {
	if [ $# -lt 6 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_SPAdes_INP_DIR=$2
				shift; shift;
				;;
			-p) NUMCPU=$2
				shift; shift;
				;;
			-m) ngsLocal_SPAdes_MAX_RAM=$2
				shift; shift;
				;;
			-k) ngsLocal_SPAdes_KMERS=$2
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
# SPAdes command, version 2.2.1
##########################################################################################

ngsCmd_SPAdes() {
	if $SE; then prnCmd "# BEGIN: SPAdes SINGLE-END ASSEMBLY"
	else prnCmd "# BEGIN: SPAdes PAIRED-END ASSEMBLY"; fi
	
	# make relevant directory
	if [ ! -d $SAMPLE/SPAdes ]; then 
		prnCmd "mkdir $SAMPLE/SPAdes"
		if ! $DEBUG; then mkdir $SAMPLE/SPAdes; fi
	fi
	
	# print version info in journal file
	prnCmd "# SPAdes v2.2.1 (check SPAdes log file)"
	
	if $SE; then
		# single-end
		prnCmd "spades.py -1 $SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_1.fq -t $NUMCPU -m $ngsLocal_SPAdes_MAX_RAM -k $ngsLocal_SPAdes_KMERS -n $SAMPLE -o $SAMPLE/SPAdes"
		if ! $DEBUG; then 
			spades.py -1 $SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_1.fq -t $NUMCPU -m $ngsLocal_SPAdes_MAX_RAM -k $ngsLocal_SPAdes_KMERS -n $SAMPLE -o $SAMPLE/SPAdes
		fi
	else
		# paired-end
		prnCmd "spades.py -1 $SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_1.fq -2 $SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_2.fq -t $NUMCPU -m $ngsLocal_SPAdes_MAX_RAM -k $ngsLocal_SPAdes_KMERS -n $SAMPLE -o $SAMPLE/SPAdes"
		if ! $DEBUG; then 
			spades.py -1 $SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_1.fq -2 $SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_2.fq -t $NUMCPU -m $ngsLocal_SPAdes_MAX_RAM -k $ngsLocal_SPAdes_KMERS -n $SAMPLE -o $SAMPLE/SPAdes
		fi
	fi
		
	# run error checking
	if ! $DEBUG; then ngsErrorChk_SPAdes $@; fi
	
	if $SE; then prnCmd "# FINISHED: SPAdes SINGLE-END ASSEMBLY"
	else prnCmd "# FINISHED: SPAdes PAIRED-END ASSEMBLY"; fi
}

##########################################################################################
# ERROR CHECKING. Make sure output files exist and are not empty.
##########################################################################################

ngsErrorChk_SPAdes() {
	prnCmd "# SPAdes ERROR CHECKING: RUNNING"

	inputFile_1="$SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_1.fq"
	inputFile_2="$SAMPLE/$ngsLocal_SPAdes_INP_DIR/unaligned_2.fq"
	outputFile="$SAMPLE/SPAdes/$SAMPLE.fasta"

	# make sure expected output file exists and is not empty
	if [ ! -s $outputFile ]; then
		errorMsg="Expected output file does not exist or is empty.\n"
		errorMsg+="\tinput file: $inputFile_1\n"
		if ! $SE; then errorMsg+="\tinput file: $inputFile_2\n"; fi
		errorMsg+="\toutput file: $outputFile\n"
		prnError "$errorMsg"
	fi

	prnCmd "# SPAdes ERROR CHECKING: DONE"
}

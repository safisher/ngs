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
# INPUT: $RAW/$SAMPLE/*_R1_*.gz
# OUTPUT: $SAMPLE/raw/unaligned_1.fq
#
# PAIRED-END READS:
# INPUT: $RAW/$SAMPLE/*_R1_*.gz and $RAW/$SAMPLE/*_R2_*.gz
# OUTPUT: $SAMPLE/raw/unaligned_1.fq and $SAMPLE/raw/unaligned_2.fq
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_INIT="Usage: `basename $0` init OPTIONS sampleID    --  prepare read file(s) for processing\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_INIT="Usage:\n\t`basename $0` init [-i inputDir] [-se] sampleID\n"
ngsHelp_INIT+="Input:\n\tinputDir/sampleID/*_R1_*.gz\n\tinputDir/sampleID/*_R2_*.gz (paired-end reads)\n"
ngsHelp_INIT+="Output:\n\tsampleID/orig/unaligned_1.fq\n\tsampleID/orig/unaligned_2.fq (paired-end reads)\n"
ngsHelp_INIT+="Options:\n"
ngsHelp_INIT+="\t-i - parent directory containing subdirectory with compressed fastq files (default: ./raw). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie inputDir/sampleID).\n"
ngsHelp_INIT+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_INIT+="By default this expects the directory './raw/sampleID' that contains the demultiplexed reads. The demultiplexed reads need to be gzipped. The files containing the first reads need to include '_R1_' in their filenames and the second read files need to contain '_R2_' in the filenames. If inputDir is used then the read files are expected to reside in 'inputDir/sampleID'. \n\n"
ngsHelp_INIT+="This will uncompress the raw files and place them in the directory './sampleID/orig'. Output files are named 'unaligned_1.fq' (first reads) and 'unaligned_2.fq' (second reads). Only unaligned_1.fq will be generated in the case of single-end reads."

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_INIT_INP_DIR=$RAW

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# INIT args: -se (optional), sampleID
##########################################################################################

ngsArgs_INIT() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_INIT_INP_DIR=$2
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
##########################################################################################

ngsCmd_INIT() {
	if $SE; then prnCmd "# BEGIN: INIT SINGLE-END"
	else prnCmd "# BEGIN: INIT PAIRED-END"; fi
	
    # make relevant directory
	if [ ! -d $SAMPLE/orig ]; then 
		prnCmd "mkdir $SAMPLE/orig"
		if ! $DEBUG; then mkdir $SAMPLE/orig; fi 
	fi
	
    # unzip raw files to orig subdirectory. Assumes contains
    # the original compressed raw files
	prnCmd "zcat $ngsLocal_INIT_INP_DIR/$SAMPLE/*_R1_* > $SAMPLE/orig/unaligned_1.fq"
	if ! $DEBUG; then 
		zcat $ngsLocal_INIT_INP_DIR/$SAMPLE/*_R1_* > $SAMPLE/orig/unaligned_1.fq
	fi
	
	if ! $SE; then 
        # paired-end
		prnCmd "zcat $ngsLocal_INIT_INP_DIR/$SAMPLE/*_R2_* > $SAMPLE/orig/unaligned_2.fq;"
		if ! $DEBUG; then 
			zcat $ngsLocal_INIT_INP_DIR/$SAMPLE/*_R2_* > $SAMPLE/orig/unaligned_2.fq
		fi
	fi
	
	if $SE; then prnCmd "# FINISHED: INIT SINGLE-END"
	else prnCmd "# FINISHED: INIT PAIRED-END"; fi
}

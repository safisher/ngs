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
# INPUT: $SRC/$SAMPLE/*R1*.gz
# OUTPUT: $SAMPLE/raw/unaligned_1.fq
#
# PAIRED-END READS:
# INPUT: $SRC/$SAMPLE/*R1*.gz and $SRC/$SAMPLE/*R2*.gz
# OUTPUT: $SAMPLE/raw/unaligned_1.fq and $SAMPLE/raw/unaligned_2.fq
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_INIT="Usage: `basename $0` init OPTIONS sampleID    --  prepare read file(s) for processing\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_INIT="Usage:\n\t`basename $0` init [-se] sampleID\n"
ngsHelp_INIT+="Input:\n\t$SRC/sampleID/*R1*.gz\n\t$SRC/sampleID/*R2*.gz (paired-end reads)\n"
ngsHelp_INIT+="Output:\n\tsampleID/orig/unaligned_1.fq\n\tsampleID/orig/unaligned_2.fq (paired-end reads)\n"
ngsHelp_INIT+="Options:\n"
ngsHelp_INIT+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_INIT+="Expects a directory called 'raw' that contains a subdirectory with demultiplexed reads. The subdirectory's name should be the sample ID. The demultiplexed reads need to be gzipped. The files containing the first reads need to include 'R1' in their filenames and the second read files need to contain 'R2' in the filenames.\n\n"
ngsHelp_INIT+="This will uncompress the raw files and place them in a subdirectory called 'orig'. Output files are named 'unaligned_1.fq' (first reads) and 'unaligned_2.fq' (second reads). Only unaligned_1.fq will be generated in the case of single-end reads."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# INIT args: -se (optional), sampleID
##########################################################################################

ngsArgs_INIT() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
	SAMPLE=$1
	
	if [ "$SAMPLE" = "-se" ]; then
	    # got -se flag instead of sample ID
		SE=true
		
	    # make sure we still have another argument, which will be the sampleID
		if [ $# -lt 2 ]; then
			printHelp $COMMAND
			exit 0
		else
			SAMPLE=$2
		fi
	fi
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
	prnCmd "zcat $SRC/$SAMPLE/*R1* > $SAMPLE/orig/unaligned_1.fq"
	if ! $DEBUG; then 
		zcat $SRC/$SAMPLE/*R1* > $SAMPLE/orig/unaligned_1.fq
	fi
	
	if ! $SE; then 
        # paired-end
		prnCmd "zcat $SRC/$SAMPLE/*R2* > $SAMPLE/orig/unaligned_2.fq;"
		if ! $DEBUG; then 
			zcat $SRC/$SAMPLE/*R2* > $SAMPLE/orig/unaligned_2.fq
		fi
	fi
	
	if $SE; then prnCmd "# FINISHED: INIT SINGLE-END"
	else prnCmd "# FINISHED: INIT PAIRED-END"; fi
}

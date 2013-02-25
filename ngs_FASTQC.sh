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
# INPUT: $SAMPLE/orig/unaligned_1.fq and $SAMPLE/trimAT/unaligned_1.fq (if present)
# OUTPUT: $SAMPLE/fastqc/* and $SAMPLE/fastqc.trim (if $SAMPLE/trimAT exists)
# REQUIRES: FastQC
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_FASTQC="Usage: `basename $0` fastqc sampleID    --  run FastQC\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_FASTQC="Usage:\n\t`basename $0` fastqc sampleID\n"
ngsHelp_FASTQC+="Input:\n\tsampleID/orig/unaligned_1.fq\n\tsampleID/trimAT/unaligned_1.fq (if sampleID/trimAT exists)\n"
ngsHelp_FASTQC+="Output:\n\tsampleID/fastqc/*\n\tsampleID/fastqc.trim (if sampleID/trimAT exists)\n"
ngsHelp_FASTQC+="Requires:\n\tFastQC ( http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ )\n\n"
ngsHelp_FASTQC+="Run FastQC on sampleID/orig/unaligned_1.fq file. If the subdirectory 'trimAT' exists (ie data has already been trimmed) then FastQC will also run on trimAT/unaligned_1.fq. When data is trimmed FastQC will automatically be run on the trimmed data, assuming the subdirectory 'fastqc' exists (ie FastQC was previously run on the untrimmed data). The FastQC output from orig will be placed in 'fastqc' while the output from trimAT will be placed in 'fastqc.trim'."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# FASTQC args: sampleID
##########################################################################################

ngsArgs_FASTQC() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	else
		SAMPLE=$1
	fi
}

##########################################################################################
# RUNNING COMMAND ACTION
# Run FASTQC on the untrimmed reads and, if present, on the trimmed reads.
##########################################################################################

ngsCmd_FASTQC() {
	prnCmd "# BEGIN: FASTQC"
	
    # make relevant directory
	if [ ! -d $SAMPLE/fastqc ]; then 
		prnCmd "mkdir $SAMPLE/fastqc"
		if ! $DEBUG; then mkdir $SAMPLE/fastqc; fi
	fi
	
    # output version of fastqc to journal
	prnCmd "# `fastqc -v`"
	
    # run fastqc on the untrimmed data
	prnCmd "fastqc --OUTDIR=$SAMPLE/fastqc $SAMPLE/orig/unaligned_1.fq"
	if ! $DEBUG; then
	    # fastqc hangs when extracting the zip file, so we do the
	    # extraction manually
		fastqc --OUTDIR=$SAMPLE/fastqc $SAMPLE/orig/unaligned_1.fq
		
	    # do some cleanup of the output files
		prnCmd "rm $SAMPLE/fastqc/unaligned_1.fq_fastqc.zip"
		rm $SAMPLE/fastqc/unaligned_1.fq_fastqc.zip
		
		prnCmd "mv $SAMPLE/fastqc/unaligned_1.fq_fastqc/* $SAMPLE/fastqc/."
		mv $SAMPLE/fastqc/unaligned_1.fq_fastqc/* $SAMPLE/fastqc/.
		
		prnCmd "rmdir $SAMPLE/fastqc/unaligned_1.fq_fastqc"
		rmdir $SAMPLE/fastqc/unaligned_1.fq_fastqc
	fi
	
    # if data has already been trimmed, then run fastqc on the
    # trimmed data
	if [ -d $SAMPLE/trimAT ]; then 
		if [ ! -d $SAMPLE/fastqc.trim ]; then 
			prnCmd "mkdir $SAMPLE/fastqc.trim"
			if ! $DEBUG; then mkdir $SAMPLE/fastqc.trim; fi
		fi
		prnCmd "fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq"
		if ! $DEBUG; then 
			fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq
			
	        # do some cleanup of the output files
			prnCmd "rm $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip"
			rm $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip
			
			prnCmd "mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/."
			mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/.
			
			prnCmd "rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc"
			rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc
		fi
	fi
	
	prnCmd "# FINISHED: FASTQC"
}

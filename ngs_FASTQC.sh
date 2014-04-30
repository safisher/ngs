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
# INPUT: $SAMPLE/init/unaligned_1.fq
# OUTPUT: $SAMPLE/fastqc.init/*
# REQUIRES: FastQC
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` fastqc OPTIONS sampleID    --  run FastQC\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_FASTQC() {
	echo -e "Usage:\n\t`basename $0` fastqc [-i inputDir] [-o outputDir] sampleID"
	echo -e "Input:\n\tsampleID/inputDir/unaligned_1.fq"
	echo -e "Output:\n\tsampleID/outputDir/*"
	echo -e "Requires:\n\tFastQC ( http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ )"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source file (default: init)."
	echo -e "\t-i outputDir - location of output files (default: fastqc). If this is changed from the default, then it will not be accessible by the STATS module.\n"
	echo -e "Run FastQC on sampleID/inputDir/unaligned_1.fq file. FastQC only uses one input file so the unaligned_1.fq file is used whether the data is single- or pair-end."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_FASTQC_INP_DIR="init"
ngsLocal_FASTQC_OUT_DIR="fastqc"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# FASTQC args: sampleID
##########################################################################################

ngsArgs_FASTQC() {
	if [ $# -lt 1 ]; then printHelp "FASTQC"; fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_FASTQC_INP_DIR=$2
				shift; shift;
				;;
			-o) ngsLocal_FASTQC_OUT_DIR=$2
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
# Run FASTQC on the untrimmed reads and, if present, on the trimmed reads.
##########################################################################################

ngsCmd_FASTQC() {
	prnCmd "# BEGIN: FASTQC"
	
    # make relevant directory
	if [ ! -d $SAMPLE/$ngsLocal_FASTQC_OUT_DIR ]; then 
		prnCmd "mkdir $SAMPLE/$ngsLocal_FASTQC_OUT_DIR"
		if ! $DEBUG; then mkdir $SAMPLE/$ngsLocal_FASTQC_OUT_DIR; fi
	fi
	
    # print version info in $SAMPLE directory
	prnCmd "# fastqc -v | awk '{print \$2}' | sed s/v//"
	if ! $DEBUG; then 
		# gets this: "-Xmx250m-Dfastqc.show_version=true-Djava.awt.headless=trueFastQC v0.10.1"
		# returns this: "0.10.1"
		ver=$(fastqc -v | awk '{print $2}' | sed s/v//)
		prnVersion "$ngsLocal_FASTQC_OUT_DIR" "program\tversion" "fastqc\t$ver"
	fi
	
    # run fastqc
	prnCmd "fastqc --OUTDIR=$SAMPLE/$ngsLocal_FASTQC_OUT_DIR $SAMPLE/$ngsLocal_FASTQC_INP_DIR/unaligned_1.fq"
	if ! $DEBUG; then
	    # fastqc hangs when extracting the zip file, so we do the
	    # extraction manually
		fastqc --OUTDIR=$SAMPLE/$ngsLocal_FASTQC_OUT_DIR $SAMPLE/$ngsLocal_FASTQC_INP_DIR/unaligned_1.fq
		
	    # do some cleanup of the output files
		prnCmd "rm $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc.zip"
		rm $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc.zip
		
		prnCmd "mv $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc/* $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/."
		mv $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc/* $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/.
		
		prnCmd "rmdir $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc"
		rmdir $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc
	fi
	
	prnCmd "# FINISHED: FASTQC"
}

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
	echo -e "Usage:\n\t`basename $0` fastqc [-i inputDir] [-o outputDir] [-se] sampleID"
	echo -e "Input:\n\tsampleID/inputDir/unaligned_1.fq"
	echo -e "Output:\n\tsampleID/outputDir/*"
	echo -e "Requires:\n\tFastQC ( http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ )"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source file (default: init)."
	echo -e "\t-i outputDir - location of output files (default: fastqc). If this is changed from the default, then it will not be accessible by the STATS module."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "Run FastQC on sampleID/inputDir/unaligned_1.fq and if PE also on sampleID/inputDir/unaligned_2.fq. FastQC only uses one input file so it is run twice in the case of pair-end reads."
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
	    # run fastqc on first reads
	    fastqc --OUTDIR=$SAMPLE/$ngsLocal_FASTQC_OUT_DIR $SAMPLE/$ngsLocal_FASTQC_INP_DIR/unaligned_1.fq
	    
	    # do some cleanup of the output files
	    # THE FOLLOWING CLEANUP IS ONLY RELEVANT TO FASTQC VERSION 0.11.1 AND LATER
	    prnCmd "mv $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc.html $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/$SAMPLE.read_1.$ngsLocal_FASTQC_OUT_DIR.html"
	    mv $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_1.fq_fastqc.html $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/$SAMPLE.read_1.$ngsLocal_FASTQC_OUT_DIR.html
	    
	    if ! $SE; then
		# run fastqc on second reads
		fastqc --OUTDIR=$SAMPLE/$ngsLocal_FASTQC_OUT_DIR $SAMPLE/$ngsLocal_FASTQC_INP_DIR/unaligned_2.fq
		
		# do some cleanup of the output files
		# THE FOLLOWING CLEANUP IS ONLY RELEVANT TO FASTQC VERSION 0.11.1 AND LATER
		prnCmd "mv $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_2.fq_fastqc.html $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/$SAMPLE.read_2.$ngsLocal_FASTQC_OUT_DIR.html"
		mv $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/unaligned_2.fq_fastqc.html $SAMPLE/$ngsLocal_FASTQC_OUT_DIR/$SAMPLE.read_2.$ngsLocal_FASTQC_OUT_DIR.html
	    fi
	fi
	
	prnCmd "# FINISHED: FASTQC"
}

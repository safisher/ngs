#!/bin/bash

# Copyright (c) 2012-2014, Stephen Fisher, Hoa Giang, and Junhyong Kim, University of
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
# INPUT: see individual commands
# OUTPUT: see individual commands
# REQUIRES: see individual commands
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` pipeline OPTIONS sampleID    --  run full pipeline\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_PIPELINE() {
	echo -e "Usage:\n\t`basename $0` pipeline [-i inputDir] [-o outputDir] [-t RNASeq | RNASeq_Stranded | RNASeq_Human | WGS] -p numProc -s species [-se] sampleID"
	echo -e "Input:\n\tsee individual commands"
	echo -e "Output:\n\tsee individual commands"
	echo -e "Requires:\n\tsee individual commands"
	echo -e "OPTIONS:"
	echo -e "\t-i - parent directory containing subdirectory with compressed fastq files (default: ./raw). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie inputDir/sampleID)."
	echo -e "\t-o - directory containing subdirectory with analysis files (default: ./analyzed). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie outputDir/sampleID)."
	echo -e "\t-t type - RNASeq or WGS (Whole Genome Sequencing) (default: RNASeq). RNASeq_Stranded assumes stranded reads for HTSeq counting and will generate intron counts. RNASeq_Human is the same as RNASeq_Stranded but also uses 'gene_name' for the name of the features in the HTSeq GTF file."
	echo -e "\t-p numProc - number of cpu to use."
	echo -e "\t-s species - species from repository: $REPO_LOCATION."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "This will process sequencing data using either an RNASeq or WGS (Whole Genome Sequencing) pipeline. For RNASeq the modules used are: init, fastqc, blast, trim, star, post, htseq, blastdb, and rsync. For WGS the modules used are: init, fastqc, blast, trim, bowtie, SPAdes, post, and rsync. See individual modules for documentation."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_PIPELINE_TYPE="RNASeq"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# PIPELINE args: -i value, -o value, -t value, -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_PIPELINE() {
	if [ $# -lt 5 ]; then printHelp "PIPELINE"; fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) RAW=$2
				shift; shift;
				;;
			-o) ANALYZED=$2
				shift; shift;
				;;
			-t) ngsLocal_PIPELINE_TYPE=$2
				shift; shift;
				;;
			-p) NUMCPU=$2
				shift; shift;
				;;
			-s) SPECIES=$2
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
# PIPELINE does not have its own command function. Rather includes the
# command functions from the following commands: init, fastqc, blast,
# trim, star, post. blastdb, htseq, rsync. See the individual
# config files.
##########################################################################################

ngsCmd_PIPELINE() {
	if $SE; then prnCmd "# BEGIN: SINGLE-END, $ngsLocal_PIPELINE_TYPE PIPELINE"
	else prnCmd "# BEGIN: PAIRED-END, $ngsLocal_PIPELINE_TYPE PIPELINE"; fi

	# PIPELINE runs the command functions from the following commands. THE
	# PIPELINE COMMAND IS SENSITIVE TO THE ORDER OF THE COMMAND FUNCTIONS
	# BELOW. For example, INIT needs to prepare the files prior to FASTQC
	# running.

	# In order to change a "default" value, the ngsArgs_XXX() function
	# needs to be called prior to the ngsCmd_XXX(). It is important
	# that $SAMPLE is included as the last argument, every time
	# ngsArgs_XXX() is called.

	########################################################
	### The modules below are included in all pipelines. ###

	# The value of $RAW is hardcoded in ngs.sh and is used to set
	# inputDir in INIT. We allow users to change this value using
	# the optional inputDir argument (-i inputDir). Since INIT
	# defaults to the original (hardcoded) value of $RAW, we need
	# to call ngsArgs_INIT() to update the value, prior to calling
	# ngsCmd_INIT().
	ngsArgs_INIT -i $RAW $SAMPLE
	ngsCmd_INIT
	ngsCmd_FASTQC
	ngsCmd_BLAST
	########################################################
	
	if [[ "$ngsLocal_PIPELINE_TYPE" = "RNASeq" ]]; then
		ngsArgs_TRIM -m 20 -q 53 -rAT 26 -rN -c $REPO_LOCATION/trim/contaminants.fa $SAMPLE
		ngsCmd_TRIM
		# Need different args to run FastQC on the trimmed data, so adjust
		# args by calling ngsArgs_FASTQC() prior to running ngsCmd_FASTQC().
		ngsArgs_FASTQC -i trim -o fastqc.trim $SAMPLE
		ngsCmd_FASTQC
		ngsCmd_STAR
		ngsCmd_HTSEQ
		#ngsCmd_BLASTDB

	elif [[ "$ngsLocal_PIPELINE_TYPE" = "RNASeq_Stranded" ]]; then
		ngsArgs_TRIM -m 20 -q 53 -rAT 26 -rN -c $REPO_LOCATION/trim/contaminants.fa $SAMPLE
		ngsCmd_TRIM
		ngsArgs_FASTQC -i trim -o fastqc.trim $SAMPLE
		ngsCmd_FASTQC
		ngsCmd_STAR
		ngsArgs_HTSEQ -stranded -introns $SAMPLE
		ngsCmd_HTSEQ

	elif [[ "$ngsLocal_PIPELINE_TYPE" = "RNASeq_Human" ]]; then
		ngsArgs_TRIM -m 20 -q 53 -rAT 26 -rN -c $REPO_LOCATION/trim/contaminants.fa $SAMPLE
		ngsCmd_TRIM
		ngsArgs_FASTQC -i trim -o fastqc.trim $SAMPLE
		ngsCmd_FASTQC
		ngsCmd_STAR
		ngsArgs_HTSEQ -stranded -introns -id gene_name $SAMPLE
		ngsCmd_HTSEQ

	elif [[ "$ngsLocal_PIPELINE_TYPE" = "WGS" ]]; then
		# disable poly-A/T trimming for WGS
		ngsArgs_TRIM -m 20 -rAT 0 -rN -c $REPO_LOCATION/trim/contaminantsMITO.fa $SAMPLE
		ngsCmd_TRIM
		ngsArgs_FASTQC -i trim -o fastqc.trim $SAMPLE
		ngsCmd_FASTQC
		ngsCmd_BOWTIE
		ngsCmd_SNP
		ngsCmd_SPADES
		ngsArgs_POST -i bowtie $SAMPLE
		ngsCmd_POST
		ngsArgs_POST -i bowtie/SE_mapping $SAMPLE
		ngsCmd_POST

	else
		prnCmd "ERROR: Invalid PIPELINE type $$ngsLocal_PIPELINE_TYPE. Valid options are 'RNASeq' and 'WGS'."
	fi

	########################################################
	### The modules below are included in all pipelines. ###

	# compress the trimmed files
	ngsCmd_POST
	# OutputDir defaults to $ANALYZED which is hardcoded in
	# ngs.sh, just like inputDir and $RAW.
	ngsArgs_RSYNC -o $ANALYZED $SAMPLE
	ngsCmd_RSYNC
	########################################################

	if $SE; then prnCmd "# FINISHED: SINGLE-END, $ngsLocal_PIPELINE_TYPE PIPELINE"
	else prnCmd "# FINISHED: PAIRED-END, $ngsLocal_PIPELINE_TYPE PIPELINE"; fi
}

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

ngsUsage_PIPELINE="Usage: `basename $0` pipeline OPTIONS sampleID    --  run full pipeline\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_PIPELINE="Usage:\n\t`basename $0` pipeline [-i inputDir] [-o outputDir] [-t RNASeq | WGS] -p numProc -s species [-se] sampleID\n"
ngsHelp_PIPELINE+="Input:\n\tsee individual commands\n"
ngsHelp_PIPELINE+="Output:\n\tsee individual commands\n"
ngsHelp_PIPELINE+="Requires:\n\tsee individual commands\n"
ngsHelp_PIPELINE+="OPTIONS:\n"
ngsHelp_PIPELINE+="\t-i - parent directory containing subdirectory with compressed fastq files (default: ./raw). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie inputDir/sampleID).\n"
ngsHelp_PIPELINE+="\t-o - directory containing subdirectory with analysis files (default: ./analyzed). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie outputDir/sampleID).\n"
ngsHelp_PIPELINE+="\t-t type - RNASeq or WGS (Whole Genome Sequencing) (default: RNASeq).\n"
ngsHelp_PIPELINE+="\t-p numProc - number of cpu to use.\n"
ngsHelp_PIPELINE+="\t-s species - species from repository: $REPO_LOCATION.\n"
ngsHelp_PIPELINE+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_PIPELINE+="This will process sequencing data using either an RNASeq or WGS (Whole Genome Sequencing) pipeline. For RNASeq the modules used are: init, fastqc, blast, trim, star, post, htseq, blastdb, and rsync. For WGS the modules used are: init, fastqc, blast, trim, bowtie, SPAdes, post, and rsync. See individual modules for documentation."

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
	if [ $# -lt 5 ]; then
		printHelp $COMMAND
		exit 0
	fi

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

	if [[ "$ngsLocal_PIPELINE_TYPE" = "RNASeq" ]]; then
		# The value of $RAW is hardcoded in ngs.sh and is used to set
		# inputDir in INIT. We allow users to change this value using the
		# optional inputDir argument (-i inputDir). Since INIT defaults to
		# the original (hardcoded) value of $RAW, we need to call
		# ngsArgs_INIT() to update the value, prior to calling
		# ngsCmd_INIT().
		ngsArgs_INIT -i $RAW $SAMPLE
		ngsCmd_INIT
		ngsCmd_FASTQC
		ngsCmd_BLAST
		ngsCmd_TRIM
		# Need different args to run FastQC on the trimmed data, so adjust
		# args by calling ngsArgs_FASTQC() prior to running ngsCmd_FASTQC().
		ngsArgs_FASTQC -i trim -o trim.fastqc $SAMPLE
		ngsCmd_FASTQC
		ngsCmd_STAR
		ngsCmd_HTSEQ
		ngsCmd_POST
		ngsCmd_BLASTDB
		# OutputDir defaults to $ANALYZED which is hardcoded in
		# ngs.sh, just like inputDir and $RAW.
		ngsArgs_RSYNC -o $ANALYZED $SAMPLE
		ngsCmd_RSYNC

	elif [[ "$ngsLocal_PIPELINE_TYPE" = "WGS" ]]; then
		ngsArgs_INIT -i $RAW $SAMPLE
		ngsCmd_INIT
		ngsCmd_FASTQC
		ngsCmd_BLAST
		ngsArgs_TRIM -c $REPO_LOCATION/trim/contaminantsMITO.fa $SAMPLE
		ngsCmd_TRIM
		ngsArgs_FASTQC -i trim -o trim.fastqc $SAMPLE
		ngsCmd_FASTQC
		ngsCmd_BOWTIE
		ngsCmd_SNP
		ngsCmd_SPAdes
		ngsCmd_POST
		ngsArgs_POST -i bowtie $SAMPLE
		ngsCmd_POST
		ngsArgs_POST -i bowtie/SE_mapping $SAMPLE
		ngsCmd_POST
		ngsArgs_RSYNC -o $ANALYZED $SAMPLE
		ngsCmd_RSYNC

	else
		prnCmd "ERROR: Invalid PIPELINE type $$ngsLocal_PIPELINE_TYPE. Valid options are 'RNASeq' and 'WGS'."
	fi

	if $SE; then prnCmd "# FINISHED: SINGLE-END, $ngsLocal_PIPELINE_TYPE PIPELINE"
	else prnCmd "# FINISHED: PAIRED-END, $ngsLocal_PIPELINE_TYPE PIPELINE"; fi
}

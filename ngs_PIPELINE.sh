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

ngsHelp_PIPELINE="Usage:\n\t`basename $0` pipeline -p numProc -s species [-se] sampleID\n"
ngsHelp_PIPELINE+="Input:\n\tsee individual commands\n"
ngsHelp_PIPELINE+="Output:\n\tsee individual commands\n"
ngsHelp_PIPELINE+="Requires:\n\tsee individual commands\n"
ngsHelp_PIPELINE+="OPTIONS:\n"
ngsHelp_PIPELINE+="\t-p numProc - number of cpu to use.\n"
ngsHelp_PIPELINE+="\t-s species - species from repository: $REPO_LOCATION.\n"
ngsHelp_PIPELINE+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_PIPELINE+="This will run init, fastqc, blast, trim, star, post, htseq, blastdb, and rsync."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# PIPELINE args: -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_PIPELINE() {
	if [ $# -lt 5 ]; then
		printHelp $COMMAND
		exit 0
	fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
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
	if $SE; then prnCmd "# BEGIN: SINGLE-END PIPELINE"
	else prnCmd "# BEGIN: PAIRED-END PIPELINE"; fi
		
	# PIPELINE runs the command functions from the following commands. THE
	# PIPELINE COMMAND IS SENSITIVE TO THE ORDER OF THE COMMAND FUNCTIONS
	# BELOW. For example, INIT needs to prepare the files prior to FASTQC
	# running.

	ngsCmd_INIT
	ngsCmd_FASTQC
	ngsCmd_BLAST
	ngsCmd_TRIM
	# need different args to run FastQC on the trimmed data, so adjust
	# args by calling ngsArgs_FASTQC() prior to running ngsCmd_FASTQC().
	ngsArgs_FASTQC -i trim -o trim.fastqc
	ngsCmd_FASTQC
	ngsCmd_STAR
	ngsCmd_POST
	ngsCmd_BLASTDB
	ngsCmd_HTSEQ
	ngsCmd_RSYNC

	if $SE; then prnCmd "# FINISHED: SINGLE-END PIPELINE"
	else prnCmd "# FINISHED: PAIRED-END PIPELINE"; fi
}

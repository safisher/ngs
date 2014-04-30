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
# INPUT: n/a
# OUTPUT: all subdirectories and files in $SAMPLE except 'init' are copied to $outputDir\$SAMPLE
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` rsync OPTIONS sampleID    --  copy data to analyzed directory\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_RSYNC() {
	echo -e "Usage:\n\t`basename $0` rsync [-o outputDir] sampleID"
	echo -e "Output\n\tall subdirectories in sampleID except 'init' are copied to $ANALYZED\sampleID"
	echo -e "Options:"
	echo -e "\t-o - directory containing subdirectory with analysis files (default: ./analyzed). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie outputDir/sampleID).\n"
	echo -e "Copies all data to 'outputDir/sampleID' directory. This will not copy the files in the sampleID/init directory. Rsync is run twice as a consistency check."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_RSYNC_OUT_DIR=$ANALYZED

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# RSYNC args: sampleID
##########################################################################################

ngsArgs_RSYNC() {
	if [ $# -lt 1 ]; then printHelp "RSYNC"; fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-o) ngsLocal_RSYNC_OUT_DIR=$2
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
# Copy results to repository.
##########################################################################################

ngsCmd_RSYNC() {
	prnCmd "# BEGIN: COPYING TO REPO"
	
    # we exclude init since that's the unaligned data which is already in the repo. Note that $JOURNAL contains the name of the log file as well as the sample directory
	prnCmd "rsync -avh --stats --exclude init --exclude $JOURNAL $SAMPLE $ngsLocal_RSYNC_OUT_DIR/."
	if ! $DEBUG; then 
		rsync -avh --stats --exclude init --exclude $JOURNAL $SAMPLE $ngsLocal_RSYNC_OUT_DIR/.
	fi
	
	# the prnCmd is here for proper annotation in the log file. The
	# command is at the end of the function so that we capture all of
	# the prnCmd output.
	prnCmd "cat $JOURNAL >> $ngsLocal_RSYNC_OUT_DIR/$JOURNAL"

	prnCmd "# FINISHED: COPYING TO REPO"
	
	# we run rsync again here to make sure everything copied
	if ! $DEBUG; then 
		rsync -avh --stats --exclude init --exclude $JOURNAL $SAMPLE $ngsLocal_RSYNC_OUT_DIR/.
	fi

	# don't squash existing journal file, so cat onto existing file
	# rather than rsync. This command is at the end of the file since
	# any prnCmd commands would not be captured.
	if ! $DEBUG; then 
		cat $JOURNAL >> $ngsLocal_RSYNC_OUT_DIR/$JOURNAL
	fi
}

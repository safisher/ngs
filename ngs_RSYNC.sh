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
# OUTPUT: all subdirectories and files in $SAMPLE except 'orig' are copied to $ANALYZED\$SAMPLE
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_RSYNC="Usage: `basename $0` rsync sampleID    --  copy data to analyzed directory\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_RSYNC="Usage:\n\t`basename $0` rsync sampleID\n"
ngsHelp_RSYNC+="Output\n\tall subdirectories in sampleID except 'orig' are copied to $ANALYZED\sampleID\n\n"
ngsHelp_RSYNC+="Copies all data to $ANALYZED directory. Does not copy sampleID/orig directory. Rsync is run twice as a consistency check."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# RSYNC args: sampleID
##########################################################################################

ngsArgs_RSYNC() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	else
		SAMPLE=$1
	fi
}

##########################################################################################
# RUNNING COMMAND ACTION
# Copy results to repository.
##########################################################################################

ngsCmd_RSYNC() {
	prnCmd "# BEGIN: COPYING TO REPO"
	
    # we exclude orig since that's the unaligned data which is already in the repo. Note that $JOURNAL contains the name of the log file as well as the sample directory
	prnCmd "rsync -avh --stats --exclude orig --exclude $JOURNAL $SAMPLE analyzed/."
	if ! $DEBUG; then 
		rsync -avh --stats --exclude orig --exclude $JOURNAL $SAMPLE analyzed/.
	fi
	
	# the prnCmd is here for proper annotation in the log file. The
	# command is at the end of the function so that we capture all of
	# the prnCmd output.
	prnCmd "cat $JOURNAL >> analyzed/$JOURNAL"

	prnCmd "# FINISHED: COPYING TO REPO"
	
	# we run rsync again here to make sure everything copied
	if ! $DEBUG; then 
		rsync -avh --stats --exclude orig --exclude $JOURNAL $SAMPLE analyzed/.
	fi

	# don't squash existing journal file, so cat onto existing file
	# rather than rsync. This command is at the end of the file since
	# any prnCmd commands would not be captured.
	if ! $DEBUG; then 
		cat $JOURNAL >> analyzed/$JOURNAL
	fi
}

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
# OUTPUT: all subdirectories except raw are copied to $ANALYZED
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_RSYNC="Usage: `basename $0` rsync sampleID    --  copy data to analyzed directory\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_RSYNC="Usage: `basename $0` rsync sampleID\n"
ngsHelp_RSYNC+="\tCopies all data to analyzed directory. Does not copy raw or bowtie directories. Rsync is run twice as a consistency check."

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
	
    # we exclude raw since that's the unaligned data which is already in the repo
	prnCmd "rsync -avh --stats --exclude raw $SAMPLE analyzed/."
	if ! $DEBUG; then 
		rsync -avh --stats --exclude raw $SAMPLE analyzed/.
	fi
	
	prnCmd "# FINISHED: COPYING TO REPO"
	
	# we run rsync again here to make sure everything copied and to
	# have the full transfer since we just changed $JOURNAL with the
	# prnCmd above.
	if ! $DEBUG; then 
		rsync -avh --stats --exclude raw $SAMPLE analyzed/.
	fi
}

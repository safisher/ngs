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
# OUTPUT: n/a
# REQUIRES: RUM version 2
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` rumstatus sampleID    --  get status of RUM run\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_RUMSTATUS() {
	echo -e "Usage:\n\t`basename $0` rumstatus sampleID\n"
	echo -e "Returns the output from RUM status command. This will only work if rumalign has previously been run."
}

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# RUMSTATUS args: sampleID
##########################################################################################

ngsArgs_RUMSTATUS() {
	if [ $# -lt 1 ]; then printHelp "RUMSTATUS"; fi

	SAMPLE=$1
}

##########################################################################################
# RUNNING COMMAND ACTION
# Return status of RUM job.
##########################################################################################

ngsCmd_RUMSTATUS() {
	prnCmd "# BEGIN: RUM STATUS"
	
	prnCmd "rum_runner status --output $SAMPLE/rum"
	if ! $DEBUG; then 
		rum_runner status --output $SAMPLE/rum
	fi
		
	prnCmd "# FINISHED: RUM STATUS"
}

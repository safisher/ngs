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
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` help COMMAND    --  expanded command help (ex: `basename $0` help blast)\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_HELP() {
	echo -e "Usage:\n\t`basename $0` help COMMAND\n"
	echo -e "Example:\n\t`basename $0` help blastdb"
}

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# HELP args: command
##########################################################################################

ngsArgs_HELP() {
	if [ $# -lt 1 ]; then printHelp "HELP"; fi
	
	# convert module name to uppercase
	helpModule=$( echo $1 | tr "[a-z]" "[A-Z]" )
	printHelp $helpModule

	exit 0
}

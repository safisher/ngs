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

ngsUsage_HELP="Usage: `basename $0` help COMMAND    --  more information (ex: `basename $0` help blast)\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_HELP="Usage: `basename $0` help COMMAND\n"
ngsHelp_HELP+="\tExample: `basename $0` help blastdb"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# HELP args: command
##########################################################################################

ngsArgs_HELP() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
	else
		printHelp $1
	fi
	exit 0
}

##########################################################################################
# RUNNING COMMAND ACTION
# There is no command action for the HELP module.
##########################################################################################


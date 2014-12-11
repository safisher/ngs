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

NGS_USAGE+="Usage: `basename $0` version    --  print version information\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_VERSION() {
	echo -e "Usage:\n\t`basename $0` version\n"
	echo -e "Prints version information for `basename $0`."
}

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# VERSION args: n/a
##########################################################################################

ngsArgs_VERSION() {
	echo "`basename $0` version $VERSION"
	exit 0
}

##########################################################################################
# RUNNING COMMAND ACTION
# There is no command action for the VERSION module.
##########################################################################################


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
# INPUT: $SAMPLE/blast/speciesCounts.txt, $SAMPLE/trim/stats.txt, $SAMPLE/star/Log.final.out
# OUTPUT: printing to console
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` stats OPTIONS modules sampleID    --  print stats from user-specified list of modules\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_STATS() {
	echo -e "Usage:\n\t`basename $0` stats [-v] modules sampleID"
	echo -e "Input:\n\tmodule specific"
	echo -e "Output:\n\tprinting to console"
	echo -e "Options:"
	echo -e "\t-v - verbose printing of alignment stats (default: off)."
	echo -e "\tmodules - comma separated, ordered list of modules to include for stats (must not include spaces).\n"
	echo -e "Prints stats for user-specified list of modules. The stats will be tab delimited so they can be copy-pasted into an Excel table. This will not write to the analysis.log file."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

# we assume this should only print the essential data as tab-delimited key:value pairs
ngsLocal_STATS_VERBOSE="false"
ngsLocal_STATS_MODULES=""

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# STATS args: sampleID
##########################################################################################

ngsArgs_STATS() {
	if [ $# -lt 2 ]; then printHelp "STATS"; fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-v) ngsLocal_STATS_VERBOSE=true
				shift;
				;;
			-*) printf "Illegal option: '%s'\n" "$1"
				printHelp $COMMAND
				exit 0
				;;
 			*) break ;;
		esac
	done

	ngsLocal_STATS_MODULES=$1

	SAMPLE=$2
}

##########################################################################################
# RUNNING COMMAND ACTION
# Print out stats
##########################################################################################

ngsCmd_STATS() {
	if $ngsLocal_STATS_VERBOSE; then
		# don't write to log file as this is just printing stats
		echo -e "\n##################################################################"
		echo -ne "# BEGIN: STATS - $ngsLocal_STATS_MODULES\n# "
		echo -ne `date`
		echo -ne "  "
		echo $SAMPLE
		echo "##################################################################"
	fi

	# these will be tab-delimited lists that include the stats column names and values
	header="Sample ID"
	values="$SAMPLE"

	# the bash IFS variable dictates the word delimiting which is "
	# \t\n" by default. We are capturing command output that includes
	# tab characters. If we don't remove the '\t' from IFS then all
	# tabs in our output get converted to spaces.
	local IFS=" "
	
	# convert module names to upper case
	ngsLocal_STATS_MODULES=$( echo $ngsLocal_STATS_MODULES | tr "[a-z]" "[A-Z]" )

	# step through each user-specified module and call the module's
	# ngsStats_MODULE() function.
	for module in ${ngsLocal_STATS_MODULES//,/ }; do
		if $ngsLocal_STATS_VERBOSE; then echo "# EXTRACTING $module STATS"; fi
		header="$header\t$(ngsStats_$module 'header')"
		values="$values\t$(ngsStats_$module 'values')"
	done

	echo -e $header
	echo -e $values

	if $ngsLocal_STATS_VERBOSE; then
		# don't write to log file as this is just printing stats
		echo -e "\n##################################################################"
		echo -ne "# FINISHED: STATS\n# "
		echo -ne `date`
		echo -ne "  "
		echo $SAMPLE
		echo "##################################################################"
	fi
}

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

ngsUsage_STATS="Usage: `basename $0` stats OPTIONS modules sampleID    --  print stats from user-specified list of modules\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_STATS="Usage:\n\t`basename $0` stats [-v] modules sampleID\n\n"
ngsHelp_STATS+="Input:\n\tsampleID/blast/speciesCounts.txt\n\tsampleID/trim/stats.txt\n\tsampleID/star/Log.final.out (if STAR specified)\n\tsampleID/rum/mapping_stats.txt (if RUM specied)\n"
ngsHelp_STATS+="Output:\n\tprinting to console\n"
ngsHelp_STATS+="Options:\n"
ngsHelp_STATS+="\t-v - verbose printing of alignment stats (default: off).\n"
ngsHelp_STATS+="\tmodules - comma separated, ordered list of modules to include for stats (must not include spaces).\n"
ngsHelp_STATS+="Prints out stats for user-specified list of modules. The stats will be tab delimited so they can be copy-pasted into an Excel table. This will not write to the analysis.log file."

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

	# we need to still have two parameters
	if [ $# -lt 2 ]; then
		printHelp $COMMAND
		exit 0
	fi

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

	for module in ${ngsLocal_STATS_MODULES//,/ }; do
		case $module in
			blast)
				if $ngsLocal_STATS_VERBOSE; then echo "# EXTRACTING BLAST STATS"; fi
				header="$header\t$(ngsStats_BLAST 'header')"
				values="$values\t$(ngsStats_BLAST 'values')"
				;;
			trim)
				if $ngsLocal_STATS_VERBOSE; then echo "# EXTRACTING TRIM STATS"; fi
				header="$header\t$(ngsStats_TRIM 'header')"
				values="$values\t$(ngsStats_TRIM 'values')"
				;;
			star)
				if $ngsLocal_STATS_VERBOSE; then echo "# EXTRACTING STAR STATS"; fi
				header="$header\t$(ngsStats_STAR 'header')"
				values="$values\t$(ngsStats_STAR 'values')"
				;;
			rum)
				if $ngsLocal_STATS_VERBOSE; then echo "# EXTRACTING RUM STATS"; fi
				header="$header\t$(ngsStats_RUM 'header')"
				values="$values\t$(ngsStats_RUM 'values')"
				;;
 			*) 
				errorMsg="Illegal module option: '$module'\n"
				prnError "$errorMsg"
				;;
		esac
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

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

ngsUsage_STATS="Usage: `basename $0` stats OPTIONS sampleID    --  print stats from blast, trimming and aligner\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_STATS="Usage:\n\t`basename $0` stats [-a aligner] sampleID\n\n"
ngsHelp_STATS+="Input:\n\tsampleID/blast/speciesCounts.txt\n\tsampleID/trim/stats.txt\n\tsampleID/star/Log.final.out\n"
ngsHelp_STATS+="Output:\n\tprinting to console\n"
ngsHelp_STATS+="Options:\n"
ngsHelp_STATS+="\t-a aligner - which aligner log file to parse for stats (default: star). Choices are 'star' or 'rum'.\n"
ngsHelp_STATS+="Prints out BLAST, TRIM, and STAR/RUM stats. The stats will be tab delimited so they can be copy-pasted into an Excel table. This will not write to the analysis.log file."

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_STATS_ALIGNER="star"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# STATS args: sampleID
##########################################################################################

ngsArgs_STATS() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-a) ngsLocal_STATS_ALIGNER=$2
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
# Print out stats
##########################################################################################

ngsCmd_STATS() {
	# don't write to log file as this is just printing stats
	echo -e "\n##################################################################"
	echo "# BEGIN: STATS"
	echo -ne `date`
	echo -ne "  "
	echo $SAMPLE
	echo "##################################################################"
	
	echo "-- BLAST --"
	cat $SAMPLE/blast/speciesCounts.txt; 

	# the last line of the speciesCounts.txt file is the following tab-delimited list of values:
	# Total Hits	Hits Not Counted	Bacteria	Fish	Fly	Human	Mouse	Rat	Yeast
	BLAST_HEADER=`tail -2 $SAMPLE/blast/speciesCounts.txt | head -1`
	BLAST_VALUES=`tail -1 $SAMPLE/blast/speciesCounts.txt`
	
	echo -e "\n-- Trimming --"
	cat $SAMPLE/trim/stats.txt; 

	TRIM_HEADER=`tail -2 $SAMPLE/trim/stats.txt | head -1`
	TRIM_VALUES=`tail -1 $SAMPLE/trim/stats.txt`

	if [ "$ngsLocal_STATS_ALIGNER" = "rum" ]; then
		echo -e "\n-- RUM --"
		#head -32 $SAMPLE/rum.trim/mapping_stats.txt

		# we aren't currently using NUMTRIM. For now it's just to flag whether or not the sample is single- or paired-end.
		numTrim=`head -32 $SAMPLE/rum.trim/mapping_stats.txt | grep "Number of read pairs:" | awk '{print $5}'`

		rumAli=`head -32 $SAMPLE/rum.trim/mapping_stats.txt | tail -10 | grep "At least one" | awk '{print $10}' | tr -d '()'`
		$RUM_HEADER="RUM Total Aligned"
		$RUM_VALUES="$rumAli"
	
		# if $NUMTRIM is empty then probably single-end, which has different stats in mapping_stats.txt file
		if [ ! "$numTrim" ]; then
			echo -e "\nPROCESSING RUM AS SINGLE-END"

			rumAli=`head -20 $SAMPLE/rum.trim/mapping_stats.txt | grep "TOTAL:" | awk '{print $3}' | tr -d '()'`
			$RUM_HEADER="RUM Total Aligned"
			$RUM_VALUES="$rumAli"
		fi

	elif [ "$ngsLocal_STATS_ALIGNER" = "star" ]; then
		echo -e "\n-- STAR --"

		avgReadLen=`grep "Average input read length" $SAMPLE/star/Log.final.out | awk -F $'\t' '{print $2}'`
		STAR_HEADER="Avg Inp Read Len"
		STAR_VALUES="$avgReadLen"

		avgMapLen=`grep "Average mapped length" $SAMPLE/star/Log.final.out | awk -F $'\t' '{print $2}'`
		STAR_HEADER="$STAR_HEADER\tAvg Uniq Map Len"
		STAR_VALUES="$STAR_VALUES\t$avgMapLen"

		uniqMap=`grep "Uniquely mapped reads %" $SAMPLE/star/Log.final.out | awk -F $'\t' '{print $2}'`
		STAR_HEADER="$STAR_HEADER\tUniq Map"
		STAR_VALUES="$STAR_VALUES\t$uniqMap"

		multimapped=`grep "% of reads mapped to multiple loci" $SAMPLE/star/Log.final.out | awk -F $'\t' '{print $2}'`
		STAR_HEADER="$STAR_HEADER\tMultimapped"
		STAR_VALUES="$STAR_VALUES\t$multimapped"

		tooShort=`grep "% of reads unmapped: too short" $SAMPLE/star/Log.final.out | awk -F $'\t' '{print $2}'`
		STAR_HEADER="$STAR_HEADER\tToo Short"
		STAR_VALUES="$STAR_VALUES\t$tooShort"
	fi

	echo -e "\n-- $SAMPLE --"

	if [ "$ngsLocal_STATS_ALIGNER" = "rum" ]; then
		echo -e "$BLAST_HEADER\t$TRIM_HEADER\t$RUM_HEADER"
		echo -e "$BLAST_VALUES\t$TRIM_VALUES\t$RUM_VALUES"

	elif [ "$ngsLocal_STATS_ALIGNER" = "star" ]; then
		echo -e "$BLAST_HEADER\t$TRIM_HEADER\t$STAR_HEADER"
		echo -e "$BLAST_VALUES\t$TRIM_VALUES\t$STAR_VALUES"

	else
		echo -e "$BLAST_HEADER\t$TRIM_HEADER"
		echo -e "$BLAST_VALUES\t$TRIM_VALUES"
	fi
	
	# don't write to log file as this is just printing stats
	echo -e "\n##################################################################"
	echo "# FINISHED: STATS"
	echo -ne `date`
	echo -ne "  "
	echo $SAMPLE
	echo "##################################################################"
}

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
# INPUT: $SAMPLE/blast/species.txt, $SAMPLE/trimAdapters.stats.txt, $SAMPLE/trimPolyAT.stats.txt, $SAMPLE/rum.trim/mappings_stats.txt
# OUTPUT: printing to console
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_STATS="Usage: `basename $0` stats sampleID    --  print stats from blast, trimming and RUM\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_STATS="Usage:\n\t`basename $0` stats sampleID\n\n"
ngsHelp_STATS+="Prints out blast, trim, and RUM stats."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# STATS args: sampleID
##########################################################################################

ngsArgs_STATS() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	else
		SAMPLE=$1
	fi
}

##########################################################################################
# RUNNING COMMAND ACTION
# Print out stats
##########################################################################################

ngsCmd_STATS() {
	prnCmd "# BEGIN: STATS"
	
	echo "-- BLAST --"
	cat $SAMPLE/blast/species.txt; 

	# the last line of the species.txt file is the following tab-delimited list of values:
	# Total Hits	Hits Not Counted	Bacteria	Fish	Fly	Human	Mouse	Rat	Yeast
	BLASTFIELDS=`tail -1 $SAMPLE/blast/species.txt`
	
	echo -e "\n-- Adapter Trimming --"
	cat $SAMPLE/trimAdapter.stats.txt; 

	NUMREADS=`grep "total reads processed" $SAMPLE/trimAdapter.stats.txt | awk '{print $1}'`
	TRIMAD_K=`grep "trimmed and kept" $SAMPLE/trimAdapter.stats.txt | awk '{print $1}'`
	TRIMAD_D=`grep "discarded with final" $SAMPLE/trimAdapter.stats.txt | awk '{print $1}'`
	
	echo -e "\n-- Poly A/T Trimming --"
	cat $SAMPLE/trimPolyAT.stats.txt; 
	
	TRIMAT_K=`grep "trimmed and kept" $SAMPLE/trimPolyAT.stats.txt | awk '{print $1}'`
	TRIMAT_D=`grep "discarded with final" $SAMPLE/trimPolyAT.stats.txt | awk '{print $1}'`
	
	echo -e "\n-- RUM --"
	head -32 $SAMPLE/rum.trim/mapping_stats.txt
	
	NUMTRIM=`head -32 $SAMPLE/rum.trim/mapping_stats.txt | grep "Number of read pairs:" | awk '{print $5}'`
	RUMALI=`head -32 $SAMPLE/rum.trim/mapping_stats.txt | tail -10 | grep "At least one" | awk '{print $10}' | tr -d '()'`
	
	# if $NUMTRIM is empty then probably single-end, which has different stats in mapping_stats.txt file
	if [ ! "$NUMTRIM" ]; then
		NUMTRIM=`head -20 $SAMPLE/rum.trim/mapping_stats.txt | grep "Number of reads:" | awk '{print $4}'`
		RUMALI=`head -20 $SAMPLE/rum.trim/mapping_stats.txt | grep "TOTAL:" | awk '{print $3}' | tr -d '()'`
		echo -e "\nPROCESSING RUM AS SINGLE-END"
	fi
	
	echo -e "\n-- $SAMPLE --"
	echo -e "Num Reads\tNum Reads Trimmed\tRUM\tTotal Hits\tHits Not Counted\tBacteria\tFish\tFly\tHuman\tMouse\tRat\tYeast\tBowtie\tAdapt Trim Kept\tAdapt Trim Disc\tPolyAT Trim Kept\tPolyAT Trim Disc\tDate"
	# extra tab between blast and trimming to account for lacking bowtie stats
	echo -e "$NUMREADS\t$NUMTRIM\t$RUMALI\t$BLASTFIELDS\t\t$TRIMAD_K\t$TRIMAD_D\t$TRIMAT_K\t$TRIMAT_D\t$(date +%m/%d/%Y)\n"
	
	prnCmd "# FINISHED: STATS"
}

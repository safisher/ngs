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
# SINGLE-END READS
# INPUT: raw/unaligned_1.fq
# OUTPUT: blastdb/*
#
# PAIRED-END READS
# INPUT: raw/unaligned_1.fq and raw/unaligned_2.fq
# OUTPUT: blastdb/*
#
# REQUIRES: makeblastdb (provided with Blast version 2)
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_BLASTDB="Usage: `basename $0` blastdb sampleID    --  create blast database from reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_BLASTDB="Usage: `basename $0` blastdb [-se] sampleID\n"
ngsHelp_BLASTDB+="\tGenerates blast database from reads.\n"
ngsHelp_BLASTDB+="\tOPTIONS:\n"
ngsHelp_BLASTDB+="\t\t-se - single-end reads (default: paired-end)"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# BLASTDB args: -se (optional), sampleID
##########################################################################################

ngsArgs_BLASTDB() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
	SAMPLE=$1
	
	if [ "$SAMPLE" = "-se" ]; then
		# got -se flag instead of sample ID
		SE=true
		
		# make sure we still have another argument, which will be the sampleID
		if [ $# -lt 2 ]; then
			printHelp $COMMAND
			exit 0
		else
			SAMPLE=$2
		fi
	fi
}

##########################################################################################
# RUNNING COMMAND ACTION
# Generate blast database from reads.
##########################################################################################

ngsCmd_BLASTDB() {
	if $SE; then prnCmd "# FINISHED: CREATE BLAST DATABASE SINGLE-END"
	else prnCmd "# FINISHED: CREATE BLAST DATABASE PAIRED-END"; fi
	
	# make relevant directory
	if [ ! -d $SAMPLE/blast.db ]; then 
		prnCmd "mkdir $SAMPLE/blast.db"
		if ! $DEBUG; then mkdir $SAMPLE/blast.db; fi
	fi
	
	# Convert raw fastq files into single fasta file
	prnCmd "awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}' $SAMPLE/raw/unaligned_1.fq > $SAMPLE/blast.db/raw.fa"
	if ! $DEBUG; then 
		awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $SAMPLE/raw/unaligned_1.fq > $SAMPLE/blast.db/raw.fa
	fi
	if ! $SE; then
		# only necessary for paired-end
		prnCmd "awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}' $SAMPLE/raw/unaligned_2.fq >> $SAMPLE/blast.db/raw.fa"
		if ! $DEBUG; then 
			awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $SAMPLE/raw/unaligned_2.fq >> $SAMPLE/blast.db/raw.fa
		fi
	fi
	
	# need to save current directory so we can return here. We also
	# need to adjust $JOURNAL so prnCmd() still works when we change
	# directories
	prnCmd "CUR_DIR=`pwd`"
	CUR_DIR=`pwd`
	
	prnCmd "JOURNAL_SAV=$JOURNAL"
	JOURNAL_SAV=$JOURNAL
	
	prnCmd "cd $SAMPLE/blast.db"
	if ! $DEBUG; then 
		cd $SAMPLE/blast.db
		
		JOURNAL=../../$JOURNAL
		prnCmd "JOURNAL=../../$JOURNAL"
	fi
	
	prnCmd "makeblastdb -in raw.fa -dbtype nucl -out $SAMPLE -title \"$SAMPLE\""
	if ! $DEBUG; then 
		makeblastdb -in raw.fa -dbtype nucl -out $SAMPLE -title "$SAMPLE"
	fi
	
	prnCmd "rm raw.fa"
	rm raw.fa
	
	# return to proper directory and restore $JOURNAL
	prnCmd "cd $CUR_DIR"
	if ! $DEBUG; then 
		cd $CUR_DIR
		
		JOURNAL=$JOURNAL_SAV
		prnCmd "JOURNAL=$JOURNAL_SAV"
	fi
	
	if $SE; then prnCmd "# FINISHED: CREATING BLAST DATABASE SINGLE-END"
	else prnCmd "# FINISHED: CREATING BLAST DATABASE PAIRED-END"; fi
}

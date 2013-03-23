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
# SINGLE-END READS:
# INPUT: $SAMPLE/orig/unaligned_1.fq
# OUTPUT: $SAMPLE/bowtie/bowtie-sorted.bam and $SAMPLE/bowtie/stats.txt
#
# PAIRED-END READS:
# INPUT: $SAMPLE/orig/unaligned_1.fq and $SAMPLE/orig/unaligned_2.fq
# OUTPUT: $SAMPLE/bowtie/bowtie-sorted.bam and $SAMPLE/bowtie/stats.txt
#
# REQUIRES: bowtie, samtools
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_BOWTIE="Usage: `basename $0` bowtie OPTIONS sampleID    --  run bowtie on untrimmed reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_BOWTIE="Usage:\n\t`basename $0` bowtie -p numProc -s species [-se] sampleID\n"
ngsHelp_BOWTIE+="Input:\n\tsampleID/orig/unaligned_1.fq\n\tsampleID/orig/unaligned_2.fq (paired-end reads)"
ngsHelp_BOWTIE+="Output:\n\tsampleID/bowtie/bowtie-sorted.bam\n\tsampleID/bowtie/stats.txt\n"
ngsHelp_BOWTIE+="Requires:\n\tbowtie ( http://bowtie-bio.sourceforge.net/index.shtml )\n\tsamtools ( http://samtools.sourceforge.net/ )\n"
ngsHelp_BOWTIE+="Options:\n"
ngsHelp_BOWTIE+="\t-p numProc - number of cpu to use\n"
ngsHelp_BOWTIE+="\t-s species - species from repository: $BOWTIE_REPO\n"
ngsHelp_BOWTIE+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_BOWTIE+="Run bowtie on the original, untrimmed data (ie sampleID/orig). Output is placed in the directory sampleID/bowtie."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# BOWTIE args: -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_BOWTIE() {
	if [ $# -lt 5 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
   	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-p) NUMCPU=$2
				shift; shift;
				;;
			-s) SPECIES=$2
				shift; shift;
				;;
			-se) SE=true
				shift;
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
# Run BOWTIE. Run BOWTIE on untrimmed data.
##########################################################################################

ngsCmd_BOWTIE() {
	if $SE; then prnCmd "# BEGIN: BOWTIE SINGLE-END ALIGNMENT"
	else prnCmd "# BEGIN: BOWTIE PAIRED-END ALIGNMENT"; fi
	
    # print version info in journal file
	prnCmd "# `bowtie --version | head -1`"
	
    # make relevant directory
	if [ ! -d $SAMPLE/bowtie ]; then 
		prnCmd "mkdir $SAMPLE/bowtie"
		if ! $DEBUG; then mkdir $SAMPLE/bowtie; fi
	fi
	
	if $SE; then 
        # single-end
		prnCmd "bowtie -t -v 3 -a -m 10 --sam -p $NUMCPU $BOWTIE_REPO/$SPECIES $SAMPLE/orig/unaligned_1.fq $SAMPLE/bowtie/output_p.sam > $SAMPLE/bowtie/stats.txt 2>&1"
		if ! $DEBUG; then 
			bowtie -t -v 3 -a -m 10 --sam -p $NUMCPU $BOWTIE_REPO/$SPECIES $SAMPLE/orig/unaligned_1.fq $SAMPLE/bowtie/output_p.sam > $SAMPLE/bowtie/stats.txt 2>&1
		fi
	else 
        # paired-end
		prnCmd "bowtie -t -v 3 -a -m 10 --sam -p $NUMCPU $BOWTIE_REPO/$SPECIES -1 $SAMPLE/orig/unaligned_1.fq -2 $SAMPLE/orig/unaligned_2.fq $SAMPLE/bowtie/output_p.sam > $SAMPLE/bowtie/stats.txt 2>&1"
		if ! $DEBUG; then 
			bowtie -t -v 3 -a -m 10 --sam -p $NUMCPU $BOWTIE_REPO/$SPECIES -1 $SAMPLE/orig/unaligned_1.fq -2 $SAMPLE/orig/unaligned_2.fq $SAMPLE/bowtie/output_p.sam > $SAMPLE/bowtie/stats.txt 2>&1
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
		cd $SAMPLE/bowtie
		
		JOURNAL=../../$JOURNAL
		prnCmd "JOURNAL=../../$JOURNAL"
	fi
	
	prnCmd "samtools view -h -b -S -o output_p.bam output_p.sam"
	if ! $DEBUG; then samtools view -h -b -S -o output_p.bam output_p.sam; fi
	prnCmd "samtools sort output_p.bam bowtie-sorted"
	if ! $DEBUG; then samtools sort output_p.bam bowtie-sorted; fi
	prnCmd "samtools index bowtie-sorted.bam"
	if ! $DEBUG; then samtools index bowtie-sorted.bam; fi
	prnCmd "rm output_p.sam output_p.bam"
	if ! $DEBUG; then rm output_p.sam output_p.bam; fi
	
    # return to proper directory and restore $JOURNAL
	prnCmd "cd $CUR_DIR"
	if ! $DEBUG; then 
		cd $CUR_DIR
		
		JOURNAL=$JOURNAL_SAV
		prnCmd "JOURNAL=$JOURNAL_SAV"
	fi
	
	if $SE; then prnCmd "# FINISHED: BOWTIE SINGLE-END ALIGNMENT"
	else prnCmd "# FINISHED: BOWTIE PAIRED-END ALIGNMENT"; fi
}

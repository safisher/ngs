#!/bin/bash

# Copyright (c) 2013, Stephen Fisher and Junhyong Kim, University of
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
# INPUT: $SAMPLE/trim/unaligned_1.fq
# OUTPUT: $SAMPLE/star/Aligned.out.sam
#
# PAIRED-END READS
# INPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq
# OUTPUT: $SAMPLE/star/Aligned.out.sam
#
# REQUIRES: STAR version 2.3.0.1 (STAR does not have a versions option so
# the version is hardcoded in this file.
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_STAR="Usage: `basename $0` star OPTIONS sampleID    --   run STAR on trimmed reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_STAR="Usage:\n\t`basename $0` star [-i inputDir] -p numProc -s species [-se] sampleID\n"
ngsHelp_STAR+="Input:\n\tsampleID/INPUTDIR/unaligned_1.fq\n\tsampleID/INPUTDIR/unaligned_2.fq (paired-end reads)\n"
ngsHelp_STAR+="Output:\n\tsampleID/star/Aligned.out.sam\n"
ngsHelp_STAR+="Requires:\n\tSTAR ( http://code.google.com/p/rna-star )\n"
ngsHelp_STAR+="Options:\n"
ngsHelp_STAR+="\t-i inputDir - location of source files (default: trim).\n"
ngsHelp_STAR+="\t-p numProc - number of cpu to use.\n"
ngsHelp_STAR+="\t-s species - species from repository: $STAR_REPO.\n"
ngsHelp_STAR+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_STAR+="Runs STAR using the trimmed files from sampleID/trim. Output is stored in sampleID/star."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# STAR args: -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_STAR() {
	if [ $# -lt 5 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
    # default value
	INPDIR="trim"

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-p) INPDIR=$2
				shift; shift;
				;;
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
# Run STAR job, assuming STAR version 2.3
##########################################################################################

ngsCmd_STAR() {
	if $SE; then prnCmd "# BEGIN: STAR SINGLE-END ALIGNMENT"
	else prnCmd "# BEGIN: STAR PAIRED-END ALIGNMENT"; fi
	
	# make relevant directory
	if [ ! -d $SAMPLE/star ]; then 
		prnCmd "mkdir $SAMPLE/star"
		if ! $DEBUG; then mkdir $SAMPLE/star; fi
	fi
	
	# print version info in journal file
	prnCmd "# assumed to be STAR v2.3.0.1"
	
	if $SE; then
		# single-end
		prnCmd "STAR --genomeDir $RUM_REPO/$SPECIES --readFilesIn $INPDIR/unaligned_1.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx"
		if ! $DEBUG; then 
			STAR --genomeDir $RUM_REPO/$SPECIES --readFilesIn $INPDIR/unaligned_1.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx
		fi
		
		prnCmd "# FINISHED: STAR SINGLE-END ALIGNMENT"
	else
		# paired-end
		prnCmd "STAR --genomeDir $RUM_REPO/$SPECIES --readFilesIn $INPDIR/unaligned_1.fq $INPDIR/unaligned_2.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx"
		if ! $DEBUG; then 
			STAR --genomeDir $RUM_REPO/$SPECIES --readFilesIn $INPDIR/unaligned_1.fq $INPDIR/unaligned_2.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx
		fi
		
		prnCmd "# FINISHED: STAR PAIRED-END ALIGNMENT"
	fi
}

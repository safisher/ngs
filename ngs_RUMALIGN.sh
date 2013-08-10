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
# INPUT: $SAMPLE/trim/unaligned_1.fq
# OUTPUT: $SAMPLE/rum.trim/RUM.bam and $SAMPLE/rum.trim/RUM_Unique.bam
#
# PAIRED-END READS
# INPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq
# OUTPUT: $SAMPLE/rum.trim/RUM.bam and $SAMPLE/rum.trim/RUM_Unique.bam
#
# REQUIRES: samtools, RUM version 2
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_RUMALIGN="Usage: `basename $0` rumalign OPTIONS sampleID    --   run RUM on trimmed reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_RUMALIGN="Usage:\n\t`basename $0` rumalign [-i inputDir] -p numProc -s species [-se] sampleID\n"
ngsHelp_RUMALIGN+="Input:\n\tsampleID/INPUTDIR/unaligned_1.fq\n\tsampleID/INPUTDIR/unaligned_2.fq (paired-end reads)\n"
ngsHelp_RUMALIGN+="Output:\n\tsampleID/rum.trim/RUM.bam (all aligned reads)\n\tsampleID/rum.trim/RUM_Unique.bam (uniquely aligned reads)\n"
ngsHelp_RUMALIGN+="Requires:\n\tRUM ( http://cbil.upenn.edu/RUM )\n"
ngsHelp_RUMALIGN+="Options:\n"
ngsHelp_RUMALIGN+="\t-i inputDir - location of source files (default: trim).\n"
ngsHelp_RUMALIGN+="\t-p numProc - number of cpu to use.\n"
ngsHelp_RUMALIGN+="\t-s species - species from repository: $RUM_REPO.\n"
ngsHelp_RUMALIGN+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_RUMALIGN+="Runs RUM using the trimmed files from sampleID/trim. Output is stored in sampleID/rum.trim directory."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# RUMALIGN args: -i value (optional), -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_RUMALIGN() {
	if [ $# -lt 5 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
    # default value
	INP_DIR="trim"

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) INP_DIR=$2
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
# Run RUM job, assuming RUM version 2.
##########################################################################################

ngsCmd_RUMALIGN() {
	if $SE; then prnCmd "# BEGIN: RUM SINGLE-END ALIGNMENT"
	else prnCmd "# BEGIN: RUM PAIRED-END ALIGNMENT"; fi
	
	# make relevant directory
	if [ ! -d $SAMPLE/rum.trim ]; then 
		prnCmd "mkdir $SAMPLE/rum.trim"
		if ! $DEBUG; then mkdir $SAMPLE/rum.trim; fi
	fi
	
	# print version info in journal file
	prnCmd "# rum_runner version"
	if ! $DEBUG; then prnCmd "# `rum_runner version`"; fi
	
	if $SE; then
		# single-end
		prnCmd "rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$INP_DIR/unaligned_1.fq"
		if ! $DEBUG; then 
			rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$INP_DIR/unaligned_1.fq
		fi
		
		prnCmd "# FINISHED: RUM SINGLE-END ALIGNMENT"
	else
		# paired-end
		prnCmd "rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$INP_DIR/unaligned_1.fq $SAMPLE/$INP_DIR/unaligned_2.fq"
		if ! $DEBUG; then 
			rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$INP_DIR/unaligned_1.fq $SAMPLE/$INP_DIR/unaligned_2.fq
		fi
		
		# run post processing to generate necessary alignment files
		rumPostProcessing $@

		prnCmd "# FINISHED: RUM PAIRED-END ALIGNMENT"
	fi
}

rumPostProcessing() {
	prnCmd "# RUM POST PROCESSING"
	
	# need to save current directory so we can return here. We also
	# need to adjust $JOURNAL so prnCmd() still works when we change
	# directories
	prnCmd "CUR_DIR=`pwd`"
	CUR_DIR=`pwd`
	
	prnCmd "JOURNAL_SAV=$JOURNAL"
	JOURNAL_SAV=$JOURNAL
	
	prnCmd "cd $SAMPLE/rum.trim"
	if ! $DEBUG; then 
		cd $SAMPLE/rum.trim
		
		JOURNAL=../../$JOURNAL
		prnCmd "JOURNAL=../../$JOURNAL"
	fi
	
	prnCmd "rm quals.fa reads.fa"
	if ! $DEBUG; then 
		rm quals.fa reads.fa
	fi
	
	# RUM version 1 generated sorted NU and Unique files. RUM version 2 does not.
	prnCmd "gzip RUM_NU RUM_Unique"
	if ! $DEBUG; then 
		gzip RUM_NU RUM_Unique
	fi
	
	prnCmd "# converting SAM output to sorted BAM file"
	prnCmd "samtools view -h -b -S -o RUM.bam RUM.sam"
	if ! $DEBUG; then 
		samtools view -h -b -S -o RUM.bam RUM.sam
	fi
	
	prnCmd "samtools sort RUM.bam RUM.sorted"
	if ! $DEBUG; then 
		samtools sort RUM.bam RUM.sorted
	fi
	
	prnCmd "samtools index RUM.sorted.bam"
	if ! $DEBUG; then 
		samtools index RUM.sorted.bam
	fi
	
	# generate BAM file containing all uniquely mapped reads. This variant will
	# remove mitochondrial genes:
	#   samtools view -H -S RUM.sam > header.sam; grep -Pv 'chrM\t' RUM.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > RUM_Unique.bam
	prnCmd "# generating RUM_Unique.bam file"
	prnCmd "samtools view -H -S RUM.sam > header.sam"
	# (1) extract all mapped reads from SAM file, (2) filter by number of mappings, (3) add header, (4) convert to BAM
	prnCmd "samtools view -S -F 0x4 RUM.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > RUM_Unique.bam"
	#prnCmd "samtools sort RUM_Unique.bam RUM_Unique.sorted"
	#prnCmd "samtools index RUM_Unique.sorted.bam"
	#prnCmd "rm header.sam RUM_Unique.bam"
	prnCmd "rm header.sam"
	if ! $DEBUG; then 
		samtools view -H -S RUM.sam > header.sam
		samtools view -S -F 0x4 RUM.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > RUM_Unique.bam
		#samtools sort RUM_Unique.bam RUM_Unique.sorted
		#samtools index RUM_Unique.sorted.bam
		#rm header.sam RUM_Unique.bam
		rm header.sam
	fi
	
	# this might be problematic if the sorting doesn't work.
	prnCmd "rm RUM.bam RUM.sam"
	if ! $DEBUG; then 
		rm RUM.bam RUM.sam
	fi
	
	# return to proper directory and restore $JOURNAL
	prnCmd "cd $CUR_DIR"
	if ! $DEBUG; then 
		cd $CUR_DIR
		
		JOURNAL=$JOURNAL_SAV
		prnCmd "JOURNAL=$JOURNAL_SAV"
	fi
}


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
# OUTPUT: $SAMPLE/star/STAR.bam and $SAMPLE/star/STAR_Unique.bam
#
# PAIRED-END READS
# INPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq
# OUTPUT: $SAMPLE/star/STAR.bam and $SAMPLE/star/STAR_Unique.bam
#
# REQUIRES: samtools, STAR version 2.3.0.1 (STAR does not have a versions option so
# the version is hardcoded in this file)
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
ngsHelp_STAR+="Output:\n\tsampleID/star/STAR.bam (all alignments)\n\tsampleID/star/STAR_Unique.bam (uniquely aligned reads)\n"
ngsHelp_STAR+="Requires:\n\tSTAR ( http://code.google.com/p/rna-star )\n"
ngsHelp_STAR+="Options:\n"
ngsHelp_STAR+="\t-i inputDir - location of source files (default: trim).\n"
ngsHelp_STAR+="\t-p numProc - number of cpu to use.\n"
ngsHelp_STAR+="\t-s species - species from repository: $STAR_REPO.\n"
ngsHelp_STAR+="\t-se - single-end reads (default: paired-end)\n\n"
ngsHelp_STAR+="Runs STAR using the trimmed files from sampleID/trim. Output is stored in sampleID/star."

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_STAR_INP_DIR="trim"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# STAR args: -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_STAR() {
	if [ $# -lt 5 ]; then
		printHelp $COMMAND
		exit 0
	fi
	
	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_STAR_INP_DIR=$2
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
	prnCmd "# STAR v2.3.0.1 (check STAR log file)"
	
	if $SE; then
		# single-end
		prnCmd "STAR --genomeDir $STAR_REPO/$SPECIES --readFilesIn $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_1.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx"
		if ! $DEBUG; then 
			STAR --genomeDir $STAR_REPO/$SPECIES --readFilesIn $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_1.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx
		fi
		
		prnCmd "# FINISHED: STAR SINGLE-END ALIGNMENT"
	else
		# paired-end
		prnCmd "STAR --genomeDir $STAR_REPO/$SPECIES --readFilesIn $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_1.fq $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_2.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx"
		if ! $DEBUG; then 
			STAR --genomeDir $STAR_REPO/$SPECIES --readFilesIn $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_1.fq $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_2.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx
		fi
		
		# run post processing to generate necessary alignment files
		starPostProcessing $@

		prnCmd "# FINISHED: STAR PAIRED-END ALIGNMENT"
	fi
}

starPostProcessing() {
	prnCmd "# STAR POST PROCESSING"
	
	# need to save current directory so we can return here. We also
	# need to adjust $JOURNAL so prnCmd() still works when we change
	# directories
	prnCmd "CUR_DIR=`pwd`"
	CUR_DIR=`pwd`
	
	prnCmd "JOURNAL_SAV=$JOURNAL"
	JOURNAL_SAV=$JOURNAL
	
	prnCmd "cd $SAMPLE/star"
	if ! $DEBUG; then 
		cd $SAMPLE/star
		
		JOURNAL=../../$JOURNAL
		prnCmd "JOURNAL=../../$JOURNAL"
	fi

	prnCmd "# compress unmapped reads files that SAM created"
	prnCmd "gzip Unmapped.out.*"
	if ! $DEBUG; then 
		gzip Unmapped.out.*
	fi
	
	prnCmd "# converting SAM output to sorted BAM file"
	prnCmd "samtools view -h -b -S -o STAR.bam Aligned.out.sam"
	if ! $DEBUG; then 
		samtools view -h -b -S -o STAR.bam Aligned.out.sam
	fi
	
	prnCmd "samtools sort STAR.bam STAR.sorted"
	if ! $DEBUG; then 
		samtools sort STAR.bam STAR.sorted
	fi
	
	prnCmd "samtools index STAR.sorted.bam"
	if ! $DEBUG; then 
		samtools index STAR.sorted.bam
	fi
	
	# generate BAM file containing all uniquely mapped reads. This variant will
	# remove mitochondrial genes:
	#   samtools view -H -S Aligned.out.sam > header.sam; grep -Pv 'chrM\t' Aligned.out.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > STAR_Unique.bam
	prnCmd "# generating STAR_Unique.bam file"
	prnCmd "samtools view -H -S Aligned.out.sam > header.sam"
	# (1) extract all mapped reads from SAM file, (2) filter by number of mappings, (3) add header, (4) convert to BAM
	prnCmd "samtools view -S -F 0x4 Aligned.out.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > STAR_Unique.bam"
	#prnCmd "samtools sort STAR_Unique.bam STAR_Unique.sorted"
	#prnCmd "samtools index STAR_Unique.sorted.bam"
	#prnCmd "rm header.sam STAR_Unique.bam"
	prnCmd "rm header.sam"
	if ! $DEBUG; then 
		samtools view -H -S Aligned.out.sam > header.sam
		samtools view -S -F 0x4 Aligned.out.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > STAR_Unique.bam
		#samtools sort STAR_Unique.bam STAR_Unique.sorted
		#samtools index STAR_Unique.sorted.bam
		#rm header.sam STAR_Unique.bam
		rm header.sam
	fi
	
	# this might be problematic if the sorting doesn't work.
	prnCmd "rm STAR.bam Aligned.out.sam"
	if ! $DEBUG; then 
		rm STAR.bam Aligned.out.sam
	fi
	
	# return to proper directory and restore $JOURNAL
	prnCmd "cd $CUR_DIR"
	if ! $DEBUG; then 
		cd $CUR_DIR
		
		JOURNAL=$JOURNAL_SAV
		prnCmd "JOURNAL=$JOURNAL_SAV"
	fi
}

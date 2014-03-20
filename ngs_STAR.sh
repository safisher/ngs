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
# OUTPUT: $SAMPLE/star/$SAMPLE.star.sorted.bam and $SAMPLE/star/$SAMPLE.star.unique.bam
#
# PAIRED-END READS
# INPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq
# OUTPUT: $SAMPLE/star/$SAMPLE.star.sorted.bam and $SAMPLE/star/$SAMPLE.star.unique.bam
#
# REQUIRES: samtools, STAR version 2.3.0.1 (STAR does not have a versions option so
# the version is hardcoded in this file)
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` star OPTIONS sampleID    --   run STAR on trimmed reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_STAR() {
	echo -e "Usage:\n\t`basename $0` star [-i inputDir] -p numProc -s species [-se] sampleID"
	echo -e "Input:\n\tsampleID/inputDir/unaligned_1.fq\n\tsampleID/inputDir/unaligned_2.fq (paired-end reads)"
	echo -e "Output:\n\tsampleID/star/sampleID.star.sorted.bam (all alignments)\n\tsampleID/star/sampleID.star.unique.bam (uniquely aligned reads)"
	echo -e "Requires:\n\tSTAR ( http://code.google.com/p/rna-star )\n\tsamtools ( http://samtools.sourceforge.net/ )"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source files (default: trim)."
	echo -e "\t-p numProc - number of cpu to use."
	echo -e "\t-s species - species from repository: $STAR_REPO."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "Runs STAR using the trimmed files from sampleID/trim. Output is stored in sampleID/star."
	echo -e "STAR options used: --genomeLoad LoadAndRemove --outReadsUnmapped Fastx"
}

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
	if [ $# -lt 5 ]; then printHelp "STAR"; fi
	
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
	
    # print version info in $SAMPLE directory
	prnCmd "# STAR version: (check STAR log file)"
	prnCmd "# samtools version: samtools 2>&1 | grep 'Version:' | awk '{print \$2}'"
	if ! $DEBUG; then 
		ver=$(samtools 2>&1 | grep 'Version:' | awk '{print $2}')
		prnVersion "star" "program\tversion\tprogram\tversion\tspecies" "star\t2.3.0e_r291\tsamtools\t$ver\t$SPECIES"
	fi

	if $SE; then
		# single-end
		prnCmd "STAR --genomeDir $STAR_REPO/$SPECIES --readFilesIn $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_1.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx"
		if ! $DEBUG; then 
			STAR --genomeDir $STAR_REPO/$SPECIES --readFilesIn $SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_1.fq --runThreadN $NUMCPU  --genomeLoad LoadAndRemove --outFileNamePrefix $SAMPLE/star/ --outReadsUnmapped Fastx
		fi
		
		# run post processing to generate necessary alignment files
		starPostProcessing $@

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
	
	prnCmd "samtools sort STAR.bam $SAMPLE.star.sorted"
	if ! $DEBUG; then 
		samtools sort STAR.bam $SAMPLE.star.sorted
	fi
	
	prnCmd "samtools index $SAMPLE.star.sorted.bam"
	if ! $DEBUG; then 
		samtools index $SAMPLE.star.sorted.bam
	fi
	
	# generate BAM file containing all uniquely mapped reads. This variant will
	# remove mitochondrial genes:
	#   samtools view -H -S Aligned.out.sam > header.sam; $GREPP -v 'chrM\t' Aligned.out.sam | $GREPP 'IH:i:1\t' | cat header.sam - | samtools view -bS - > STAR_Unique.bam
	prnCmd "# generating STAR_Unique.bam file"
	prnCmd "samtools view -H -S Aligned.out.sam > header.sam"
	# (1) extract all mapped reads from SAM file, (2) filter by number of mappings, (3) add header, (4) convert to BAM
	prnCmd "samtools view -S -F 0x4 Aligned.out.sam | $GREPP 'NH:i:1\t' | cat header.sam - | samtools view -bS - > $SAMPLE.star.unique.bam"
	prnCmd "rm header.sam"
	if ! $DEBUG; then 
		samtools view -H -S Aligned.out.sam > header.sam
		samtools view -S -F 0x4 Aligned.out.sam | $GREPP 'NH:i:1\t' | cat header.sam - | samtools view -bS - > $SAMPLE.star.unique.bam
		rm header.sam
	fi
	
	# this might be problematic if the sorting doesn't work.
	prnCmd "rm STAR.bam Aligned.out.sam"
	if ! $DEBUG; then 
		rm STAR.bam Aligned.out.sam
	fi
	
	# rename output stats file to conform to other modules
	prnCmd "mv $SAMPLE/star/Log.final.out $SAMPLE/star/$SAMPLE.star.stats.txt"
	if ! $DEBUG; then
		mv $SAMPLE/star/Log.final.out $SAMPLE/star/$SAMPLE.star.stats.txt
	fi
	
	# return to proper directory and restore $JOURNAL
	prnCmd "cd $CUR_DIR"
	if ! $DEBUG; then 
		cd $CUR_DIR
		
		JOURNAL=$JOURNAL_SAV
		prnCmd "JOURNAL=$JOURNAL_SAV"
	fi

	# run error checking
	if ! $DEBUG; then ngsErrorChk_STAR $@; fi
}

##########################################################################################
# ERROR CHECKING. Make sure output files exist and are not empty.
##########################################################################################

ngsErrorChk_STAR() {
	prnCmd "# STAR ERROR CHECKING: RUNNING"

	inputFile_1="$SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_1.fq"
	inputFile_2="$SAMPLE/$ngsLocal_STAR_INP_DIR/unaligned_2.fq"
	outputFile_1="$SAMPLE/star/$SAMPLE.star.sorted.bam"
	outputFile_2="$SAMPLE/star/$SAMPLE.star.unique.bam"

	# make sure expected output files exists
	if [[ ! -s $outputFile_1 || ! -s $outputFile_2 ]]; then
		errorMsg="Error with STAR output files (don't exist or are empty).\n"
		errorMsg+="\tinput file: $inputFile_1\n"
		if ! $SE; then errorMsg+="\tinput file: $inputFile_2\n"; fi
		errorMsg+="\toutput file (sorted alignments): $outputFile_1\n"
		errorMsg+="\toutput file (unique alignments): $outputFile_2\n"
		prnError "$errorMsg"
	fi

	prnCmd "# STAR ERROR CHECKING: DONE"
}

##########################################################################################
# PRINT STATS. Prints a tab-delimited list stats of interest.
##########################################################################################

ngsStats_STAR() {
	if [ $# -ne 1 ]; then
		prnError "Incorrect number of parameters for ngsStats_STAR()."
	fi

	avgReadLen=`grep "Average input read length" $SAMPLE/star/$SAMPLE.star.stats.txt | awk -F $'\t' '{print $2}'`
	STAR_HEADER="Avg Inp Read Len"
	STAR_VALUES="$avgReadLen"

	avgMapLen=`grep "Average mapped length" $SAMPLE/star/$SAMPLE.star.stats.txt | awk -F $'\t' '{print $2}'`
	STAR_HEADER="$STAR_HEADER\tAvg Uniq Map Len"
	STAR_VALUES="$STAR_VALUES\t$avgMapLen"

	uniqMap=`grep "Uniquely mapped reads %" $SAMPLE/star/$SAMPLE.star.stats.txt | awk -F $'\t' '{print $2}'`
	STAR_HEADER="$STAR_HEADER\tUniq Map"
	STAR_VALUES="$STAR_VALUES\t$uniqMap"

	multimapped=`grep "% of reads mapped to multiple loci" $SAMPLE/star/$SAMPLE.star.stats.txt | awk -F $'\t' '{print $2}'`
	STAR_HEADER="$STAR_HEADER\tMultimapped"
	STAR_VALUES="$STAR_VALUES\t$multimapped"

	tooShort=`grep "% of reads unmapped: too short" $SAMPLE/star/$SAMPLE.star.stats.txt | awk -F $'\t' '{print $2}'`
	STAR_HEADER="$STAR_HEADER\tToo Short"
	STAR_VALUES="$STAR_VALUES\t$tooShort"

	case $1 in
		header)
			# the second to the last line of the stats.txt file is a tab-delimited lists of headers
			echo "$STAR_HEADER"
			;;

		values)
			# the last line of the stats.txt file is a tab-delimited lists of values
			echo "$STAR_VALUES"
			;;

		*) 
			# incorrect argument
			prnError "Invalid parameter for ngsStats_STAR() (got $1, expected: 'header|values')."
			;;
	esac
}

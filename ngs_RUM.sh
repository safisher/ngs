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
# OUTPUT: $SAMPLE/rum/$SAMPLE.rum.sorted.bam and $SAMPLE/rum/$SAMPLE.rum.unique.bam
#
# PAIRED-END READS
# INPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq
# OUTPUT: $SAMPLE/rum/$SAMPLE.rum.sorted.bam and $SAMPLE/rum/$SAMPLE.rum.unique.bam
#
# REQUIRES: samtools, RUM version 2
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` rum OPTIONS sampleID    --   run RUM on trimmed reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_RUM() {
	echo -e "Usage:\n\t`basename $0` rum [-i inputDir] -p numProc -s species [-se] sampleID"
	echo -e "Input:\n\tsampleID/INPUTDIR/unaligned_1.fq\n\tsampleID/INPUTDIR/unaligned_2.fq (paired-end reads)"
	echo -e "Output:\n\tsampleID/rum/sampleID.rum.sorted.bam (all aligned reads)\n\tsampleID/rum/sampleID.rum.unique.bam (uniquely aligned reads)"
	echo -e "Requires:\n\tRUM ( http://cbil.upenn.edu/RUM )\n\tsamtools ( http://samtools.sourceforge.net/ )"
	echo -e "Options:"
	echo -e "\t-i inputDir - location of source files (default: trim)."
	echo -e "\t-p numProc - number of cpu to use."
	echo -e "\t-s species - species from repository: $RUM_REPO."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "Runs RUM using the trimmed files from sampleID/trim. Output is stored in sampleID/rum directory. No non-default options are specified for RUM."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_RUM_INP_DIR="trim"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# RUM args: -i value (optional), -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_RUM() {
	if [ $# -lt 5 ]; then printHelp "RUM"; fi
	
	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_RUM_INP_DIR=$2
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

ngsCmd_RUM() {
	if $SE; then prnCmd "# BEGIN: RUM SINGLE-END ALIGNMENT"
	else prnCmd "# BEGIN: RUM PAIRED-END ALIGNMENT"; fi
	
	# make relevant directory
	if [ ! -d $SAMPLE/rum ]; then 
		prnCmd "mkdir $SAMPLE/rum"
		if ! $DEBUG; then mkdir $SAMPLE/rum; fi
	fi

    # print version info in $SAMPLE directory
	prnCmd "# rum_runner version | awk '{print $3}' | sed s/v// | sed s/,//"
	if ! $DEBUG; then 
		# gets this: "RUM version v2.0.3_04, released November 12, 2012"
		# returns this: "2.0.3_04"
		ver=$(rum_runner version | awk '{print $3}' | sed s/v// | sed s/,//)
		prnVersion "rum" "rum_version\tspecies" "$ver\t$SPECIES"
	fi
	
	if $SE; then
		# single-end
		prnCmd "rum_runner align --output $SAMPLE/rum --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_1.fq"
		if ! $DEBUG; then 
			rum_runner align --output $SAMPLE/rum --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_1.fq
		fi
		
		prnCmd "# FINISHED: RUM SINGLE-END ALIGNMENT"
	else
		# paired-end
		prnCmd "rum_runner align --output $SAMPLE/rum --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_1.fq $SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_2.fq"
		if ! $DEBUG; then 
			rum_runner align --output $SAMPLE/rum --name $SAMPLE --index $RUM_REPO/$SPECIES --chunks $NUMCPU $SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_1.fq $SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_2.fq
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
	
	prnCmd "cd $SAMPLE/rum"
	if ! $DEBUG; then 
		cd $SAMPLE/rum
		
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
	
	prnCmd "samtools sort RUM.bam $SAMPLE.rum.sorted"
	if ! $DEBUG; then 
		samtools sort RUM.bam $SAMPLE.rum.sorted
	fi
	
	prnCmd "samtools index $SAMPLE.rum.sorted.bam"
	if ! $DEBUG; then 
		samtools index $SAMPLE.rum.sorted.bam
	fi
	
	# generate BAM file containing all uniquely mapped reads. This variant will
	# remove mitochondrial genes:
	#   samtools view -H -S RUM.sam > header.sam; grep -Pv 'chrM\t' RUM.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > RUM_Unique.bam
	prnCmd "# generating $SAMPLE.rum.unique.bam file"
	prnCmd "samtools view -H -S RUM.sam > header.sam"
	# (1) extract all mapped reads from SAM file, (2) filter by number of mappings, (3) add header, (4) convert to BAM
	prnCmd "samtools view -S -F 0x4 RUM.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > $SAMPLE.rum.unique.bam"
	#prnCmd "samtools sort RUM_Unique.bam RUM_Unique.sorted"
	#prnCmd "samtools index RUM_Unique.sorted.bam"
	#prnCmd "rm header.sam RUM_Unique.bam"
	prnCmd "rm header.sam"
	if ! $DEBUG; then 
		samtools view -H -S RUM.sam > header.sam
		samtools view -S -F 0x4 RUM.sam | grep -P 'IH:i:1\t' | cat header.sam - | samtools view -bS - > $SAMPLE.rum.unique.bam
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

##########################################################################################
# ERROR CHECKING. Make sure output files exist and are not empty.
##########################################################################################

ngsErrorChqk_RUM() {
	prnCmd "# RUM ERROR CHECKING: RUNNING"

	inputFile_1="$SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_1.fq"
	inputFile_2="$SAMPLE/$ngsLocal_RUM_INP_DIR/unaligned_2.fq"
	outputFile_1="$SAMPLE/rum/$SAMPLE.rum.sorted.bam"
	outputFile_2="$SAMPLE/rum/$SAMPLE.rum.unique.bam"

	# make sure expected output files exists
	if [[ ! -s $outputFile_1 || ! -s $outputFile_2 ]]; then
		errorMsg="Error with RUM output files (don't exist or are empty).\n"
		errorMsg+="\tinput file: $inputFile_1\n"
		if ! $SE; then errorMsg+="\tinput file: $inputFile_2\n"; fi
		errorMsg+="\toutput file (sorted alignments): $outputFile_1\n"
		errorMsg+="\toutput file (unique alignments): $outputFile_2\n"
		prnError "$errorMsg"
	fi

	prnCmd "# RUM ERROR CHECKING: DONE"
}

##########################################################################################
# PRINT STATS. Prints a tab-delimited list stats of interest.
##########################################################################################

ngsStats_RUM() {
	if [ $# -ne 1 ]; then
		prnError "Incorrect number of parameters for ngsStats_RUM()."
	fi

	# we aren't currently using NUMTRIM. For now it's just to flag whether or not the sample is SE or PE.
	numTrim=`head -32 $SAMPLE/rum/mapping_stats.txt | grep "Number of read pairs:" | awk '{print $5}'`
	
	rumAli=`head -32 $SAMPLE/rum/mapping_stats.txt | tail -10 | grep "At least one" | awk '{print $10}' | tr -d '()'`
	# if $numTrim is empty then probably SE, which has different stats in mapping_stats.txt file
	if [ ! "$numTrim" ]; then
		rumAli=`head -20 $SAMPLE/rum/mapping_stats.txt | grep "TOTAL:" | awk '{print $3}' | tr -d '()'`
	fi

	RUM_HEADER="RUM Total Aligned"
	RUM_VALUES="$rumAli"

	case $1 in
		header)
			# the second to the last line of the stats.txt file is a tab-delimited lists of headers
			echo "$RUM_HEADER"
			;;

		values)
			# the last line of the stats.txt file is a tab-delimited lists of values
			echo "$RUM_VALUES"
			;;

		*) 
			# incorrect argument
			prnError "Invalid parameter for ngsStats_RUM() (got $1, expected: 'header|values')."
			;;
	esac
}

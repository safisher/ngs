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
# INPUT: $SAMPLE/rum.trim/RUM.sam
# OUTPUT: $SAMPLE/rum.trim/RUM.sorted.bam, $SAMPLE/rum.trim/RUM_Unique.sorted.bam
#         directory $SAMPLE/trimAT renamed $SAMPLE/trim.Ad.PolyAT and contents compressed
#         directory $SAMPLE/trimAD deleted
# REQUIRES: samtools
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

ngsUsage_POST="Usage: `basename $0` post sampleID    --  clean up RUM and trimmed data\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_POST="Usage:\n\t`basename $0` post sampleID\n"
ngsHelp_POST+="Input:\n\tsampleID/rum.trim/RUM.sam\n"
ngsHelp_POST+="Output:\n\tsampleID/rum.trim/RUM.sorted.bam\n\tsampleID/rum.trim/RUM_Unique.sorted.bam\n\tdirectory sampleID/trimAT renamed sampleID/trim.Ad.PolyAT and contents compressed\n\tdirectory sampleID/trimAD deleted\n"
ngsHelp_POST+="Requires:\n\tsamtools ( http://samtools.sourceforge.net/ )\n\n"
ngsHelp_POST+="Cleans up RUM output, compressing files as feasible, converting SAM output to sorted BAM, running parseFeatureQuant.py script to summarize transcript output, removed trimAD, compresses trimAT files and moves them into a directory called 'trim.Ad.PolyAT'."

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# POST args: sampleID
##########################################################################################

ngsArgs_POST() {
	if [ $# -lt 1 ]; then
		printHelp $COMMAND
		exit 0
	else
		SAMPLE=$1
	fi
}

##########################################################################################
# RUNNING COMMAND ACTION
# POST command. Post-processing of RUM results. Removes excess RUM and
# trim files, converts SAM to BAM, and compresses trimming data that's
# being saved.
##########################################################################################

ngsCmd_POST() {
	prnCmd "# BEGIN: POST PROCESSING"
	
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
	prnCmd "samtools view -S -F 0x4 RUM.sam | grep -P 'IH:i:1\t' RUM.sam | cat header.sam - | samtools view -bS - > RUM_Unique.bam"
	prnCmd "samtools sort RUM_Unique.bam RUM_Unique.sorted"
	prnCmd "samtools index RUM_Unique.sorted.bam"
	prnCmd "rm header.sam RUM_Unique.bam"
	if ! $DEBUG; then 
		samtools view -H -S RUM.sam > header.sam
		samtools view -S -F 0x4 RUM.sam | grep -P 'IH:i:1\t' RUM.sam | cat header.sam - | samtools view -bS - > RUM_Unique.bam
		samtools sort RUM_Unique.bam RUM_Unique.sorted
		samtools index RUM_Unique.sorted.bam
		rm header.sam RUM_Unique.bam
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
	
	prnCmd "rm $SAMPLE/trimAD/*"
	if ! $DEBUG; then 
		rm $SAMPLE/trimAD/*
	fi
	
	prnCmd "rmdir $SAMPLE/trimAD"
	if ! $DEBUG; then 
		rmdir $SAMPLE/trimAD
	fi
	
	prnCmd "gzip $SAMPLE/trimAT/*fq"
	if ! $DEBUG; then 
		gzip $SAMPLE/trimAT/*fq
	fi
	
	prnCmd "mv $SAMPLE/trimAT $SAMPLE/trim.Ad.PolyAT"
	if ! $DEBUG; then 
		mv $SAMPLE/trimAT $SAMPLE/trim.Ad.PolyAT
	fi
	
	prnCmd "# FINISHED: POST PROCESSING"
}

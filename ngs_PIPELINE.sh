#!/bin/bash

# Copyright (c) 2012-2014, Stephen Fisher, Hoa Giang, and Junhyong Kim, University of
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
# INPUT: see individual commands
# OUTPUT: see individual commands
# REQUIRES: see individual commands
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` pipeline OPTIONS sampleID    --  run full pipeline\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_PIPELINE() {
	echo -e "Usage:\n\t`basename $0` pipeline [-i inputDir] [-o outputDir] [-t RNASeq | RNASeq_Stranded | RNASeq_Human | WGS] [-l readLength] -p numProc -s species [-se] sampleID"
	echo -e "Input:\n\tsee individual commands"
	echo -e "Output:\n\tsee individual commands"
	echo -e "Requires:\n\tsee individual commands"
	echo -e "OPTIONS:"
	echo -e "\t-i - parent directory containing subdirectory with compressed fastq files (default: ./raw). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie inputDir/sampleID)."
	echo -e "\t-o - directory containing subdirectory with analysis files (default: ./analyzed). This is the parent directory of the sample-specific directory. The sampleID will be used to complete the directory path (ie outputDir/sampleID)."
	echo -e "\t-t type - RNASeq, RNASeq_BC, or WGS (Whole Genome Sequencing) (default: RNASeq). The RNASeq_BC adds barcode selection and trimming to the standard RNASeq pipeline."
	echo -e "\t-l readLength - read length (default = 100). If paired end then this is the length of one mate. Used to determine blast e-values and star library length."
	echo -e "\t-p numProc - number of cpu to use."
	echo -e "\t-s species - species from repository: $REPO_LOCATION."
	echo -e "\t-se - single-end reads (default: paired-end)\n"
	echo -e "\t-c contaminant - name of contaminant file from $REPO_LOCATION/trim to be used for trimming. Default is contaminants.fa"
	echo -e "\t-f features - list of feature types for quantification. Features will be assigned hierarchically, in the order listed. Availible features are listed in the header of GTF files at $REPO_LOCATION/verse. Default is exon."
	echo -e "\t-id id_attr - ID attribute from a GTF file you will use for quantification. Final gene counts will be output using this field. Should be either gene_name or gene_id. Default is gene_id."
	echo -e "\t-stranded - data comes from stranded library preparation. Reads will only be counted if they align on the same strand as annotated features. Default is unstranded."
	echo -e "\t-lines_sines - quantify LINE and SINE elements, separately from other features. "
	echo -e "\t-chgrp group - change the unix group of a sample and it's data before syncing to the repo. Default is no change."
	echo -e "\t-small - abreviated pipeline that keeps fewer large output files. The trim output files, the position sorted STAR output, and the STAR unmapped reads files are all excluded.\n"
	echo -e "This will process sequencing data using either an RNASeq or WGS (Whole Genome Sequencing) pipeline. For RNASeq the modules used are: init, fastqc, blast, trim, star, verse, post, and rsync. For WGS the modules used are: init, fastqc, blast, trim, bowtie, SPAdes, post, and rsync. See individual modules for documentation."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_PIPELINE_TYPE="RNASeq"
ngsLocal_CONTAM_NAME="contaminants.fa"
ngsLocal_PIPELINE_FEATURES="exon"
ngsLocal_PIPELINE_LINES_SINES=""
ngsLocal_PIPELINE_ID_ATTR="gene_id"
ngsLocal_PIPELINE_STRANDED=""
ngsLocal_PIPELINE_GROUP=""
small=false
noinit=false

#########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# PIPELINE args: -i value, -o value, -t value, -p value, -s value, -se (optional), sampleID
##########################################################################################

ngsArgs_PIPELINE() {
	if [ $# -lt 5 ]; then printHelp "PIPELINE"; fi

	ngsLocal_PIPELINE_INITIAL_ARGS="$@" #needed to recurse for multi-barcode samples

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) RAW=$2		#global
				shift; shift;
				;;
			-o) ANALYZED=$2		#global
				shift; shift;
				;;
			-l) READ_LENGTH=$2	#global
				shift; shift;
				;;
			-p) NUMCPU=$2		#global
				shift; shift;
				;;
			-s) SPECIES=$2		#global
				shift; shift;
				;;
			-se) SE=true		#global
				shift;
				;;
			-t) ngsLocal_PIPELINE_TYPE=$2
				shift; shift;
				;;
			-c) ngsLocal_CONTAM_NAME=$2
				shift; shift;
				;;
			-f) ngsLocal_PIPELINE_FEATURES=$2
				shift; shift;
				;;
			-b) ngsLocal_PIPELINE_BC=$2
				shift; shift;
				;;
			-id) ngsLocal_PIPELINE_ID_ATTR=$2
				shift; shift;
				;;
			-stranded) ngsLocal_PIPELINE_STRANDED="-stranded"
				shift;
				;;
			-lines_sines) ngsLocal_PIPELINE_LINES_SINES="-lines_sines"
				shift;
				;;
			-chgrp) ngsLocal_PIPELINE_GROUP="-g $2"
				shift; shift;
				;;
			-small) small=true
				shift;
				;;
			-noinit) noinit=true
				shift;
				;;
			-*) printf "Illegal option: '%s'\n" "$1"
				printHelp $COMMAND
				exit 0
				;;
			*) break ;;
		esac
	done
	
	SAMPLE=$1   #global
}

##########################################################################################
# RUNNING COMMAND ACTION
# PIPELINE does not have its own command function. Rather includes the
# command functions from the following commands: init, fastqc, blast,
# trim, star, post. blastdb, htseq, rsync. See the individual
# config files.
##########################################################################################

ngsCmd_PIPELINE() {
	if $SE; then prnCmd "# BEGIN: SINGLE-END, $ngsLocal_PIPELINE_TYPE PIPELINE"
	else prnCmd "# BEGIN: PAIRED-END, $ngsLocal_PIPELINE_TYPE PIPELINE"; fi

	# PIPELINE runs the command functions from the following commands. THE
	# PIPELINE COMMAND IS SENSITIVE TO THE ORDER OF THE COMMAND FUNCTIONS
	# BELOW. For example, INIT needs to prepare the files prior to FASTQC
	# running.

	# In order to change a "default" value, the ngsArgs_XXX() function
	# needs to be called prior to the ngsCmd_XXX(). It is important
	# that $SAMPLE is included as the last argument, every time
	# ngsArgs_XXX() is called.

	########################################################
	### The modules below are included in all pipelines. ###

	# The value of $RAW is hardcoded in ngs.sh and is used to set
	# inputDir in INIT. We allow users to change this value using
	# the optional inputDir argument (-i inputDir). Since INIT
	# defaults to the original (hardcoded) value of $RAW, we need
	# to call ngsArgs_INIT() to update the value, prior to calling
	# ngsCmd_INIT().

	if ! $noinit; then
	    ngsArgs_INIT -i $RAW $SAMPLE
	    ngsCmd_INIT
	fi
	ngsCmd_FASTQC
   ########################################################

      if [[ "$ngsLocal_PIPELINE_TYPE" = "RNASeq" ]]; then
	    ngsArgs_BLAST -l $READ_LENGTH -k TATAGTGAGT -p $NUMCPU -s $SPECIES $SAMPLE
	    ngsCmd_BLAST
	    ngsArgs_TRIM -t $NUMCPU -m 20 -q 53 -rAT 26 -rN -c $ngsLocal_CONTAM_NAME $SAMPLE
	    ngsCmd_TRIM
	    # Need different args to run FastQC on the trimmed data, so adjust
	    # args by calling ngsArgs_FASTQC() prior to running ngsCmd_FASTQC().
	    ngsArgs_FASTQC -i trim -o fastqc.trim $SAMPLE
	    ngsCmd_FASTQC
	    ngsCmd_STAR
	    ngsArgs_VERSE $ngsLocal_PIPELINE_STRANDED $ngsLocal_PIPELINE_LINES_SINES -l $ngsLocal_PIPELINE_FEATURES -id $ngsLocal_PIPELINE_ID_ATTR -p $NUMCPU -s $SPECIES $SAMPLE
	    ngsCmd_VERSE
	    ngsLocal_PIPELINE_check_readCnts
	    #ngsCmd_BLASTDB
                          
	elif [[ "$ngsLocal_PIPELINE_TYPE" = "RNASeq_BC" ]]; then
	    ngsArgs_BLAST -l $READ_LENGTH -k TATAGTGAGT -p $NUMCPU -s $SPECIES $SAMPLE
	    ngsCmd_BLAST
	    ngsArgs_BARCODE -b $ngsLocal_PIPELINE_BC $SAMPLE
	    ngsCmd_BARCODE
	    ngsArgs_TRIM -i barcode.trim -t $NUMCPU -rPoly -rN -m 10 -c $ngsLocal_CONTAM_NAME $SAMPLE
	    ngsCmd_TRIM
	    ngsArgs_FASTQC -i barcode.trim -o fastqc.barcode $SAMPLE
	    ngsCmd_FASTQC
	    ngsCmd_STAR
	    ngsArgs_VERSE $ngsLocal_PIPELINE_STRANDED $ngsLocal_PIPELINE_LINES_SINES -l $ngsLocal_PIPELINE_FEATURES -id $ngsLocal_PIPELINE_ID_ATTR -p $NUMCPU -s $SPECIES $SAMPLE
	    ngsCmd_VERSE
	    ngsLocal_PIPELINE_check_readCnts
	    ngsArgs_POST -i barcode.trim $SAMPLE
	    ngsCmd_POST
	    ngsArgs_POST -i trim $SAMPLE  #Args has to be called to reset prior call

	elif [[ "$ngsLocal_PIPELINE_TYPE" = "WGS" ]]; then
	    ngsArgs_BLAST -l $READ_LENGTH -k TATAGTGAGT -p $NUMCPU -s $SPECIES $SAMPLE
	    ngsCmd_BLAST
	    # disable poly-A/T trimming for WGS
	    if [[ $ngsLocal_CONTAM_NAME == "contaminants.fa" ]]; then
		# If the contaminants file is still the default, change it to the WGS default. 
		ngsLocal_CONTAM_NAME="contaminantsWGS.fa"
	    fi
	    ngsArgs_TRIM -m 19 -rAT 0 -rN -c contaminantsWGS.fa $SAMPLE
	    ngsCmd_TRIM
	    ngsArgs_FASTQC -i trim -o fastqc.trim $SAMPLE
	    ngsCmd_FASTQC
	    ngsCmd_BOWTIE
	    ngsCmd_SNP
	    ngsCmd_SPADES
	    ngsArgs_POST -i bowtie $SAMPLE
	    ngsCmd_POST
	    ngsArgs_POST -i bowtie/SE_mapping $SAMPLE
	    ngsCmd_POST
	    #ngsArgs_POST -i trim $SAMPLE  #Args has to be called to reset prior call

	else
	    prnCmd "ERROR: Invalid PIPELINE type $ngsLocal_PIPELINE_TYPE. Valid options are 'RNASeq' and 'WGS'."
	fi

	########################################################
	### The modules below are included in all pipelines. ###

	# compress the trimmed files
	if $small; then
	    prnCmd "chgrp -R $ngsLocal_PIPELINE_GROUP $SAMPLE"
	    chgrp -R repo-admin $SAMPLE

	    prnCmd "rsync -avh --stats --exclude "init" --exclude "unaligned_*.fq" --exclude "Unmapped.out.*" $SAMPLE $ANALYZED/"
	    rsync -avh --stats --exclude "init" --exclude "unaligned_*.fq" --exclude "Unmapped.out.*" $SAMPLE $ANALYZED/
	else
	    ngsArgs_POST -i trim $ngsLocal_PIPELINE_GROUP $SAMPLE #Change group to repo-admin before syncing. Uses -i trim to reset any previous calls to args
	    ngsCmd_POST
	    # OutputDir defaults to $ANALYZED which is hardcoded in
	    # ngs.sh, just like inputDir and $RAW.
	    ngsArgs_RSYNC -o $ANALYZED $SAMPLE
	    ngsCmd_RSYNC
	fi
	########################################################

	if $SE; then prnCmd "# FINISHED: SINGLE-END, $ngsLocal_PIPELINE_TYPE PIPELINE"
	else prnCmd "# FINISHED: PAIRED-END, $ngsLocal_PIPELINE_TYPE PIPELINE"; fi
}

ngsLocal_PIPELINE_check_readCnts() {
    local starPairs=$(grep "Uniquely mapped reads number" $SAMPLE/star/$SAMPLE.star.stats.txt | cut -f2)
    local versePairs=$(grep "TotalRead" $SAMPLE/verse/$SAMPLE.verse.summary.txt | cut -f2)

    if [[ $starPairs -ne $versePairs ]]; then
	prnError "Star reports $starPairs unique mapped pairs, but Verse reports $versePairs input pairs"
    fi

}

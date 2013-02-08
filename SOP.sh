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

VERSION=1.2.8
# VERSION 1.2.8, 1/27/13
#   - added running of HTSeq (runHTSeq.py)
# VERSION 1.2.7, 1/25/13
#   - updated generating of RUM_Unique so that fixmate is run and the file is sorted.
# VERSION 1.2.6, 1/12/13
#   - added creating of blast database
# VERSION 1.2.5, 12/18/12
#   - removed generate of SAMPLE.txt ouptut file since it doesn't work for non-mouse species
#   - added generate of RUM_unique.bam output file, containing all unique reads
# VERSION 1.2.4, 10/6/12
#   - tweaked stats to work with single-end
#   - added # before comment lines in sop.txt output
#   - added TODO's
# VERSION 1.2.3, 10/4/12
#   - added allprep1 and all1 to facilitate single-end processing
# VERSION 1.2.2, 9/24/12
#   - added trimming and RUM processing of single-end reads
# VERSION 1.2.1, 9/18/12
#   - updated arguments to be action specific
#   - changed input directory from unaligned to raw and output directory from aligned to analyzed
#   - no longer hard codes location of RUM binaries.
# VERSION 1.2, 9/11/12
#   - updated for RUM 2.0.2_03
#     - remove rumpre, rumali, rumpost and added rumalign and rumresume
#     - RUM option: limit-bowtie-nu is now a RUM default
# VERSION 1.1, 8/29/12
#   - added DEBUG option
#   - changed how and when subdirectories are created
#   - added more runtime options (allrum, allprep) and changed the definition of 'all'
#   - clearified console and log output
#   - removed RUM v1 code
#   - using RUM flag --limit-bowtie-nu
#   - added fastqc

# TODO:
#   - allow for options in any order
#   - create single-end flag instead of separate options

cmdOpts() {
	echo -e "\nUSAGE: SOP.sh SAMPLE ACTION OPTIONS\n"
	echo -e "ACTIONS:"
	echo -e ""
	echo -e "\tinit - uncompress raw files placing them in raw directory. Assumes compressed first read file contains R1 label and second read contains R2 label. Output files created are labeled unaligned_1.fq and unaligned_2.fq."
	echo -e ""
	echo -e "\tfastqc - run FastQC on raw/unaligned_1.fq file. If the subdirectory 'trimAT' exists (ie data has already been trimmed) then FastQC will also run on trimAT/unaligned_1.fq. When data is trimmed FastQC will automatically be run on the trimmed data, assuming the subdirectory 'fastqc' exists (ie FastQC was previously run on the untrimmed data). The FastQC output from raw will be placed in 'fastqc' while the output from trimAT will be placed in 'fastqc.trim'."
	echo -e ""
	echo -e "\tblast - run blast on 5000 reads randomly sampled from raw/unaligned_1.fq. Blast paramters used are '-num_descriptions 10 -num_alignments 10 -word_size 15 -gapopen 3 -gapextend 1 -culling_limit 1 -evalue 1e-15'. The output is put in a directory called blast."
	echo -e "\t\tOPTIONS:"
	echo -e "\t\t\tnumProc - number of cpu to use."
	echo -e ""
	echo -e "\tbowtie - run bowtie on raw data. See the script runBowtiePaired.sh for bowtie parameters used. Output is placed in the directory called bowtie. Bowtie uses raw, untrimmed, data as input."
	echo -e "\t\tOPTIONS: numProc must preceed species"
	echo -e "\t\t\tnumProc - number of cpu to use."
	echo -e "\t\t\tspecies - species files 'hg19', 'mm9', 'rat.m4' are located in /lab/repo/resources/indexes."
	echo -e ""
	echo -e "\ttrim | trim1 - runs trimAdapters.py followed by trimPolyAT.py to trim data. Adapter trimmed data is placed in trimAD while PolyAT trimmed data is placed in trimAT. Trimming must be done in order for RUM to work as it uses the files in trimAT. For single-end reads, use 'trim1' instead of 'trim'. Single-end reads are processed with trimAdaptersSingle.py and trimPolyATSingle.py."
	echo -e ""
	echo -e "\trumalign | rumalign1 - runs RUM using the trimmed files from trimAT. Output is stored in directory 'rum.trim'. For single-end reads, use 'rumalign1' instead of 'rumalign'."
	echo -e "\t\tOPTIONS: numProc must preceed species"
	echo -e "\t\t\tnumProc - number of cpu to use."
	echo -e "\t\t\tspecies - species files 'drosophila', 'hg19', 'mm9', 'rat', 'saccer3', and 'zebrafish' are located in /lab/repo/resources/rum2."
	echo -e ""
	echo -e "\trumresume - continues a halted RUM run. This will only work if rumalign has previously been run."
	echo -e ""
	echo -e "\trumstatus - returns the output from RUM status. This will only work if rumalign has previously been run."
	echo -e ""
	echo -e "\tpost - cleans up RUM output, compressing files as feasible, converting SAM output to sorted BAM, running parseFeatureQuant.py script to summarize transcript output, removed trimAD, compresses trimAT files and moves them into a directory called 'trim.Ad.PolyAT'."
	echo -e ""
	echo -e "\thtseq - run HTSeq using runHTSeq.py script. This requires the sorted BAM file containing unique reads that is generated by 'post'."
	echo -e "\t\t\tspecies - species specific library with gene model ('zebrafish')."
	echo -e "\t\t\tprefix - identifier to extract all gene IDs from output. For example 'ENSDARG' is prefix for all zebrafish genes."
	echo -e ""
	echo -e "\trsync - Copies all data to repo analyzed directory. Does not copy raw or bowtie directories. Rsync is run twice as a consistency check."
	echo -e ""
	echo -e "\tstats - prints out blast, trim, and RUM stats."
	echo -e ""
	echo -e "\tblastdb - generates blast database from reads."
	echo -e ""
	echo -e "\tallprep - includes init, fastqc, blast, and trim."
	echo -e ""
	echo -e "\tall - includes allprep, rumalign, post, htseq, and rsync."
	echo -e "\t\t\targuments - numProc (for blast and RUM), species (for RUM), library (for htSeq), prefix (for HTSeq)."
	echo -e ""
	echo -e "\tallprep1 - includes init, fastqc, blast, and trim1 (ie single-end version of allprep)."
	echo -e ""
	echo -e "\tall1 - includes allprep1, rumalign1, post, and rsync (ie single-end version of all)."
	echo -e ""
	echo -e "\thelp - this output.\n"
	echo -e ""
	echo -e "SOP.sh requires a subdirectory (or symbolic link to a directory) called 'raw'. The raw directory needs to contain a subdirectory for each sample, with each subdirectory containing the unaligned data. The first read should contain an 'R1' in the name while the second read should contain an 'R2' in the name.\n"
	echo -e "SOP.sh also requires a subdirectory (or symbolic link to a directory) called 'analyzed'. The analyzed directory is where SOP.sh will copy the aligned data.\n"
	echo -e ""
	echo "To run an initialization of a SAMPLE called 'Test' use:"
	echo "    SOP.sh Test init"
	echo -e "\nTo run RUM on a SAMPLE called 'Test' use:"
	echo "    SOP.sh Test rumalign mm9 32"
}

# make comparisons case insensitive
shopt -s nocasematch

if [[ -z "$2" || $2 == "help" ]]; then
	cmdOpts
	exit 0 
fi

SAMPLE=$1
ACTION=$2

if [[ "$ACTION" == *"all"* || "$ACTION" == "blast" || "$ACTION" == "bowtie" || "$ACTION" == "rumalign" || "$ACTION" == "rumalign1" ]]; then
	if [[ -z "$3" ]]; then
		cmdOpts
		exit 0 
	fi
	NPROC=$3
fi
if [[ "$ACTION" == "all1" || "$ACTION" == "bowtie" || "$ACTION" == "rumalign" || "$ACTION" == "rumalign1" ]]; then
	if [[ -z "$4" ]]; then
		cmdOpts
		exit 0 
	fi
	SPECIES=$4
fi
if [[ "$ACTION" == "all" ]]; then
	if [[ -z "$4" ]]; then
		cmdOpts
		exit 0 
	fi
	SPECIES=$4
fi
#if [[ "$ACTION" == "all" ]]; then
#	if [[ -z "$5" ]]; then
#		cmdOpts
#		exit 0 
#	fi
#	SPECIES=$4
#	PREFIX=$5
#fi
if [[ "$ACTION" == "htseq"  ]]; then
	if [[ -z "$4" ]]; then
		cmdOpts
		exit 0 
	fi
	SPECIES=$3
	PREFIX=$4
fi

OUT="$SAMPLE/sop.txt"

DEBUG=false   # disable commands

prnCmd() {
	if [[ $1 == *"# BEGIN:"* ]]; then
		# copy to console
		echo
		echo "##################################################################"
		echo $1
		echo -ne `date`
		echo -ne "  "
		echo $SAMPLE
		echo "##################################################################"

		echo -ne `date` >> $OUT
		echo -ne "\t" >> $OUT
		echo "##################################################################" >> $OUT
	fi
	echo -ne `date` >> $OUT
	echo -ne "\t" >> $OUT
	echo -n $1 >> $OUT
	echo >> $OUT
	if [[ $1 == *"# FINISHED"* ]]; then
		echo -ne `date` >> $OUT
		echo -ne "\t" >> $OUT
		echo "##################################################################" >> $OUT

		# insert extra line between sections
		echo >> $OUT

		# copy to console
		echo
		echo "##################################################################"
		echo $1
		echo -ne `date`
		echo -ne "  "
		echo $SAMPLE
		echo "##################################################################"
	fi
}

# this needs to happen first, so we can create the output file that
# is used by prnCmd ($OUT). This happens even during debugging because this
# directory is where the $OUT file is located by default.
if [ ! -d $SAMPLE ]; then
	mkdir $SAMPLE
fi
prnCmd "# SOP.sh Version: "$VERSION

#################################################
# Initial config
#################################################
if [[ "$ACTION" == *"all"* || "$ACTION" == "init" ]]; then
	prnCmd "# BEGIN: INIT"

	# make relevant directory
	if [ ! -d $SAMPLE/raw ]; then 
		prnCmd "mkdir $SAMPLE/raw"
		if ! $DEBUG; then mkdir $SAMPLE/raw; fi 
	fi

	# unzip raw files to raw directory. Assumes raw is a link to the directory containing the original compressed raw files
	prnCmd "zcat raw/$SAMPLE/*R1* > $SAMPLE/raw/unaligned_1.fq; zcat raw/$SAMPLE/*R2* > $SAMPLE/raw/unaligned_2.fq;"
	if ! $DEBUG; then 
		zcat raw/$SAMPLE/*R1* > $SAMPLE/raw/unaligned_1.fq; zcat raw/$SAMPLE/*R2* > $SAMPLE/raw/unaligned_2.fq
	fi

	prnCmd "# FINISHED: INIT"
fi

#################################################
# Run FASTQC. Run FASTQC on untrimmed data.
#################################################
if [[ "$ACTION" == *"all"* || "$ACTION" == "fastqc" ]]; then
	prnCmd "# BEGIN: FASTQC"

	# make relevant directory
	if [ ! -d $SAMPLE/fastqc ]; then 
		prnCmd "mkdir $SAMPLE/fastqc"
		if ! $DEBUG; then mkdir $SAMPLE/fastqc; fi
	fi

	# run fastqc on the untrimmed data
	prnCmd "fastqc --OUTDIR=$SAMPLE/fastqc $SAMPLE/raw/unaligned_1.fq"
	if ! $DEBUG; then
		# fastqc hangs when extracting the zip file, so we do the extraction manually
		fastqc --OUTDIR=$SAMPLE/fastqc $SAMPLE/raw/unaligned_1.fq

		#prnCmd "mv $SAMPLE/fastqc/unaligned_1.fq_fastqc.zip $SAMPLE/fastqc/$SAMPLE.fastqc.zip"
		#mv $SAMPLE/fastqc/unaligned_1.fq_fastqc.zip $SAMPLE/fastqc/$SAMPLE.fastqc.zip

		# do some cleanup of the output files
		prnCmd "rm $SAMPLE/fastqc/unaligned_1.fq_fastqc.zip"
		rm $SAMPLE/fastqc/unaligned_1.fq_fastqc.zip

		prnCmd "mv $SAMPLE/fastqc/unaligned_1.fq_fastqc/* $SAMPLE/fastqc/."
		mv $SAMPLE/fastqc/unaligned_1.fq_fastqc/* $SAMPLE/fastqc/.

		prnCmd "rmdir $SAMPLE/fastqc/unaligned_1.fq_fastqc"
		rmdir $SAMPLE/fastqc/unaligned_1.fq_fastqc
	fi

	# if data has already been trimmed, then run fastqc on the trimmed data
	if [ -d $SAMPLE/trimAT ]; then 
		if [ ! -d $SAMPLE/fastqc.trim ]; then 
			prnCmd "mkdir $SAMPLE/fastqc.trim"
			if ! $DEBUG; then mkdir $SAMPLE/fastqc.trim; fi
		fi
		prnCmd "fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq"
		if ! $DEBUG; then 
			fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq
			
		    # do some cleanup of the output files
			prnCmd "rm $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip"
			rm $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip
			
			prnCmd "mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/."
			mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/.
			
			prnCmd "rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc"
			rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc
		fi
	fi

	prnCmd "# FINISHED: FASTQC"
fi

#################################################
# Run BLAST. This will do a BLAST search on 5,000 
# untrimmed reads, using the nt database.
#################################################
if [[ "$ACTION" == *"all"* || "$ACTION" == "blast" ]]; then
	prnCmd "# BEGIN: BLAST"

	# make relevant directory
	if [ ! -d $SAMPLE/blast ]; then 
		prnCmd "mkdir $SAMPLE/blast"
		if ! $DEBUG; then mkdir $SAMPLE/blast; fi
	fi

	# Get 5,000 randomly sampled reads
	# Usage: randomSample.py <num lines> <lines grouped> <input> <output>
	prnCmd "randomSample.py 5000 4 $SAMPLE/raw/unaligned_1.fq $SAMPLE/blast/raw.fq > $SAMPLE/blast/species.txt"
	if ! $DEBUG; then 
		randomSample.py 5000 4 $SAMPLE/raw/unaligned_1.fq $SAMPLE/blast/raw.fq > $SAMPLE/blast/species.txt
	fi

	# Convert fastq file to fasta file
	prnCmd "awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}' $SAMPLE/blast/raw.fq > $SAMPLE/blast/raw.fa"
	if ! $DEBUG; then 
		awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $SAMPLE/blast/raw.fq > $SAMPLE/blast/raw.fa
	fi

	# Run BLAST. Output file should end with ".txt"
	prnCmd "blastn -query $SAMPLE/blast/raw.fa -db nt -num_descriptions 10 -num_alignments 10 -word_size 15 -gapopen 3 -gapextend 1 -culling_limit 1 -evalue 1e-15 -num_threads $NPROC -out $SAMPLE/blast/blast.txt"
	if ! $DEBUG; then 
		blastn -query $SAMPLE/blast/raw.fa -db nt -num_descriptions 10 -num_alignments 10 -word_size 15 -gapopen 3 -gapextend 1 -culling_limit 1 -evalue 1e-15 -num_threads $NPROC -out $SAMPLE/blast/blast.txt
	fi

	# Parse BLAST output. Adds ".txt" to input file. Will generate *.cvs and *.hits files.
	# Usage: parseBlast.py <reads fasta file> <blast prefix>
	prnCmd "parseBlast.py $SAMPLE/blast/raw.fa $SAMPLE/blast/blast >> $SAMPLE/blast/species.txt"
	if ! $DEBUG; then 
		parseBlast.py $SAMPLE/blast/raw.fa $SAMPLE/blast/blast >> $SAMPLE/blast/species.txt
	fi

	prnCmd "# FINISHED: BLAST"
fi

#################################################
# Run BOWTIE. Run BOWTIE on untrimmed data.
#################################################
if [[ "$ACTION" == "bowtie" ]]; then
	prnCmd "# BEGIN: BOWTIE"

	# make relevant directory
	if [ ! -d $SAMPLE/bowtie ]; then 
		prnCmd "mkdir $SAMPLE/bowtie"
		if ! $DEBUG; then mkdir $SAMPLE/bowtie; fi
	fi

	# runBowtiePaired.sh <input dir> <output dir> <ref seq -> hg19 | mm9 | rat.m4> <num nodes>
	prnCmd "runBowtiePaired.sh $SAMPLE/raw $SAMPLE/bowtie $SPECIES $NPROC"
	if ! $DEBUG; then 
		runBowtiePaired.sh $SAMPLE/raw $SAMPLE/bowtie $SPECIES $NPROC
	fi

	prnCmd "# FINISHED: BOWTIE"
fi

#################################################
# Trim single-end reads
#################################################
if [[ "$ACTION" == "allprep1" || "$ACTION" == "all1" || "$ACTION" == "trim1" ]]; then
	prnCmd "# BEGIN: SINGLE-END TRIMMING"

	# make relevant directory
	if ! $DEBUG; then 
		prnCmd "mkdir $SAMPLE/trimAD $SAMPLE/trimAT"
		if [ ! -d $SAMPLE/trimAD ]; then mkdir $SAMPLE/trimAD; fi
		if [ ! -d $SAMPLE/trimAT ]; then mkdir $SAMPLE/trimAT; fi
	fi
	prnCmd "trimAdaptersSingle.py $SAMPLE/raw/unaligned_1.fq $SAMPLE/trimAD > $SAMPLE/trimAdapter.stats.txt"
	if ! $DEBUG; then 
		trimAdaptersSingle.py $SAMPLE/raw/unaligned_1.fq $SAMPLE/trimAD > $SAMPLE/trimAdapter.stats.txt
	fi

	prnCmd "trimPolyATSingle.py $SAMPLE/trimAD/unaligned_1.fq $SAMPLE/trimAT > $SAMPLE/trimPolyAT.stats.txt"
	if ! $DEBUG; then 
		trimPolyATSingle.py $SAMPLE/trimAD/unaligned_1.fq $SAMPLE/trimAT > $SAMPLE/trimPolyAT.stats.txt
	fi

	# if we ran fastqc on raw data (ie $SAMPLE/fastqc exists), then
	# also run on trimmed data. Put this fastqc output in separate
	# directory so it doesn't squash the output from raw
	if [ -d $SAMPLE/fastqc ]; then 
		if [ ! -d $SAMPLE/fastqc.trim ]; then 
			prnCmd "mkdir $SAMPLE/fastqc.trim"
			if ! $DEBUG; then mkdir $SAMPLE/fastqc.trim; fi
		fi
		prnCmd "fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq"
		if ! $DEBUG; then 
			fastqc --noextract --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq

		    # do some cleanup of the output files
			prnCmd "mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip $SAMPLE/fastqc.trim/$SAMPLE.fastqc.zip"
			mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip $SAMPLE/fastqc.trim/$SAMPLE.fastqc.zip
		fi
	fi

	prnCmd "# FINISHED: SINGLE-END TRIMMING"
fi

#################################################
# Trim reads
#################################################
if [[ "$ACTION" == "allprep" || "$ACTION" == "all" || "$ACTION" == "trim" ]]; then
	prnCmd "# BEGIN: PAIRED-END TRIMMING"

	# make relevant directory
	if ! $DEBUG; then 
		prnCmd "mkdir $SAMPLE/trimAD $SAMPLE/trimAT"
		if [ ! -d $SAMPLE/trimAD ]; then mkdir $SAMPLE/trimAD; fi
		if [ ! -d $SAMPLE/trimAT ]; then mkdir $SAMPLE/trimAT; fi
	fi

	prnCmd "trimAdapters.py $SAMPLE/raw/unaligned_1.fq $SAMPLE/raw/unaligned_2.fq $SAMPLE/trimAD > $SAMPLE/trimAdapter.stats.txt"
	if ! $DEBUG; then 
		trimAdapters.py $SAMPLE/raw/unaligned_1.fq $SAMPLE/raw/unaligned_2.fq $SAMPLE/trimAD > $SAMPLE/trimAdapter.stats.txt
	fi

	prnCmd "trimPolyAT.py $SAMPLE/trimAD/unaligned_1.fq $SAMPLE/trimAD/unaligned_2.fq $SAMPLE/trimAT > $SAMPLE/trimPolyAT.stats.txt"
	if ! $DEBUG; then 
		trimPolyAT.py $SAMPLE/trimAD/unaligned_1.fq $SAMPLE/trimAD/unaligned_2.fq $SAMPLE/trimAT > $SAMPLE/trimPolyAT.stats.txt
	fi

	# if we ran fastqc on raw data (ie $SAMPLE/fastqc exists), then
	# also run on trimmed data. Put this fastqc output in separate
	# directory so it doesn't squash the output from raw
	if [ -d $SAMPLE/fastqc ]; then 
		if [ ! -d $SAMPLE/fastqc.trim ]; then 
			prnCmd "mkdir $SAMPLE/fastqc.trim"
			if ! $DEBUG; then mkdir $SAMPLE/fastqc.trim; fi
		fi
		prnCmd "fastqc --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq"
		if ! $DEBUG; then 
			fastqc --noextract --OUTDIR=$SAMPLE/fastqc.trim $SAMPLE/trimAT/unaligned_1.fq

		    # do some cleanup of the output files
			prnCmd "mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip $SAMPLE/fastqc.trim/$SAMPLE.fastqc.zip"
			mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc.zip $SAMPLE/fastqc.trim/$SAMPLE.fastqc.zip
			
			#prnCmd "mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/."
			#mv $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc/* $SAMPLE/fastqc.trim/.
			
			#prnCmd "rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc"
			#rmdir $SAMPLE/fastqc.trim/unaligned_1.fq_fastqc
		fi
	fi

	prnCmd "# FINISHED: PAIRED-END TRIMMING"
fi

#################################################
# Run RUM version 2
#################################################
if [[ "$ACTION" == "all" || "$ACTION" == "rumalign" ]]; then
	prnCmd "# BEGIN: RUM PAIRED-END ALIGNMENT"

	# make relevant directory
	if [ ! -d $SAMPLE/rum.trim ]; then 
		prnCmd "mkdir $SAMPLE/rum.trim"
		if ! $DEBUG; then mkdir $SAMPLE/rum.trim; fi
	fi

	# SPECIES is assumed to be mm9, rat, hg19, saccer3 (yeast)
	prnCmd "rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index /lab/repo/resources/rum2/$SPECIES --chunks $NPROC $SAMPLE/trimAT/unaligned_1.fq $SAMPLE/trimAT/unaligned_2.fq"
	if ! $DEBUG; then 
		rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index /lab/repo/resources/rum2/$SPECIES --chunks $NPROC $SAMPLE/trimAT/unaligned_1.fq $SAMPLE/trimAT/unaligned_2.fq
	fi

	prnCmd "# FINISHED: RUM PAIRED-END ALIGNMENT"
fi

if [[ "$ACTION" == "all1" || "$ACTION" == "rumalign1" ]]; then
	prnCmd "# BEGIN: RUM SINGLE-END ALIGNMENT"

	# make relevant directory
	if [ ! -d $SAMPLE/rum.trim ]; then 
		prnCmd "mkdir $SAMPLE/rum.trim"
		if ! $DEBUG; then mkdir $SAMPLE/rum.trim; fi
	fi

	# SPECIES is assumed to be mm9, rat, hg19, saccer3 (yeast)
	prnCmd "rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index /lab/repo/resources/rum2/$SPECIES --chunks $NPROC $SAMPLE/trimAT/unaligned_1.fq"
	if ! $DEBUG; then 
		rum_runner align --output $SAMPLE/rum.trim --name $SAMPLE --index /lab/repo/resources/rum2/$SPECIES --chunks $NPROC $SAMPLE/trimAT/unaligned_1.fq
	fi

	prnCmd "# FINISHED: RUM SINGLE-END ALIGNMENT"
fi

if [[ "$ACTION" == "rumresume" ]]; then
	prnCmd "# BEGIN: RUM RESUME"

	prnCmd "rum_runner resume --output $SAMPLE/rum.trim"
	if ! $DEBUG; then 
		rum_runner resume --output $SAMPLE/rum.trim
	fi

	prnCmd "# FINISHED: RUM RESUME"
fi

if [[ "$ACTION" == "rumstatus" ]]; then
	prnCmd "# BEGIN: RUM STATUS"

	prnCmd "rum_runner status --output $SAMPLE/rum.trim"
	if ! $DEBUG; then 
		rum_runner status --output $SAMPLE/rum.trim
	fi

	prnCmd "# FINISHED: RUM STATUS"
fi

#################################################
# Post-processing of results. Removes excess RUM 
# and trim files, converts SAM to BAM, and compresses
# trimming data that's being saved.
#################################################
if [[ "$ACTION" == "all" || "$ACTION" == "all1" || "$ACTION" == "post" ]]; then
	prnCmd "# BEGIN: POST PROCESSING"

	# need to save current directory so we can return here. We also
	# need to adjust $OUT so prnCmd() still works when we change
	# directories
	prnCmd "CUR_DIR=`pwd`"
	CUR_DIR=`pwd`

	prnCmd "OUT_SAV=$OUT"
	OUT_SAV=$OUT

	prnCmd "cd $SAMPLE/rum.trim"
	if ! $DEBUG; then 
		cd $SAMPLE/rum.trim

		OUT=../../$OUT
		prnCmd "OUT=../../$OUT"
	fi

	# this isn't working for non-mouse species.
	#prnCmd "parseFeatureQuant.py feature_quantifications_$SAMPLE $SAMPLE.txt > $SAMPLE.log.txt"
	#if ! $DEBUG; then 
		#parseFeatureQuant.py feature_quantifications_$SAMPLE $SAMPLE.txt > $SAMPLE.log.txt
	#fi

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
	prnCmd "grep -P 'IH:i:1\t' RUM.sam | cat header.sam - | samtools view -bS - > RUM_Unique.bam"
	prnCmd "samtools sort RUM_Unique.bam RUM_Unique.sorted"
	prnCmd "samtools index RUM_Unique.sorted.bam"
	prnCmd "rm header.sam RUM_Unique.bam"
	if ! $DEBUG; then 
		samtools view -H -S RUM.sam > header.sam
		grep -P 'IH:i:1\t' RUM.sam | cat header.sam - | samtools view -bS - > RUM_Unique.bam
		samtools sort RUM_Unique.bam RUM_Unique.sorted
		samtools index RUM_Unique.sorted.bam
		rm header.sam RUM_Unique.bam
	fi

	# this might be problematic if the sorting doesn't work.
	prnCmd "rm RUM.bam RUM.sam"
	if ! $DEBUG; then 
		rm RUM.bam RUM.sam
	fi

	# return to proper directory and restore $OUT
	prnCmd "cd $CUR_DIR"
	if ! $DEBUG; then 
		cd $CUR_DIR

		OUT=$OUT_SAV
		prnCmd "OUT=$OUT_SAV"
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
fi

#################################################
# Run HTSeq
#################################################
#if [[ "$ACTION" == "all" || "$ACTION" == "all1" || "$ACTION" == "htseq" ]]; then
if [[ "$ACTION" == "htseq" ]]; then
	prnCmd "# BEGIN: RUNNING HTSEQ"

	# make relevant directory
	if [ ! -d $SAMPLE/htseq ]; then 
		prnCmd "mkdir $SAMPLE/htseq"
		if ! $DEBUG; then mkdir $SAMPLE/htseq; fi
	fi

	# We assume that RUM worked and 'post' has completed.
	prnCmd "runHTSeq.py $SAMPLE/rum.trim/RUM_Unique.sorted.bam $SAMPLE/htseq/$SAMPLE /lab/repo/resources/htseq/$SPECIES/$SPECIES.gz"
	if ! $DEBUG; then 
		runHTSeq.py $SAMPLE/rum.trim/RUM_Unique.sorted.bam $SAMPLE/htseq/$SAMPLE /lab/repo/resources/htseq/$SPECIES/$SPECIES.gz
	fi

	# parse output into three files: gene counts ($SAMPLE.htseq.cnts.txt), 
	# warnings ($SAMPLE.htseq.err.txt), log ($SAMPLE.htseq.log.txt)
	prnCmd "grep $PREFIX $SAMPLE/htseq/$SAMPLE.htseq.out > $SAMPLE/htseq/$SAMPLE.htseq.cnts.txt"
	prnCmd "grep -v $PREFIX $SAMPLE/htseq/$SAMPLE.htseq.out | grep -v Warning > $SAMPLE/htseq/$SAMPLE.htseq.log.txt"
	prnCmd "grep -v $PREFIX $SAMPLE/htseq/$SAMPLE.htseq.out | grep Warning > $SAMPLE/htseq/$SAMPLE.htseq.err.txt"
	prnCmd "rm $SAMPLE/htseq/$SAMPLE.htseq.out"
	if ! $DEBUG; then 
		grep $PREFIX $SAMPLE/htseq/$SAMPLE.htseq.out > $SAMPLE/htseq/$SAMPLE.htseq.cnts.txt
		grep -v $PREFIX $SAMPLE/htseq/$SAMPLE.htseq.out | grep -v Warning > $SAMPLE/htseq/$SAMPLE.htseq.log.txt
		grep -v $PREFIX $SAMPLE/htseq/$SAMPLE.htseq.out | grep Warning > $SAMPLE/htseq/$SAMPLE.htseq.err.txt
		rm $SAMPLE/htseq/$SAMPLE.htseq.out
	fi

	prnCmd "# FINISHED: RUNNING HTSEQ"
fi

#################################################
# Copy results to repo
#################################################
if [[ "$ACTION" == "all" || "$ACTION" == "all1" || "$ACTION" == "rsync" ]]; then
	prnCmd "# BEGIN: COPYING TO REPO"

	# We assume that RUM worked and hence we don't need to store
	# bowtie output in the repo. We are including blast to give more
	# info about the reads.
	prnCmd "rsync -avh --stats --exclude raw --exclude bowtie $SAMPLE analyzed/."
	if ! $DEBUG; then 
		rsync -avh --stats --exclude raw --exclude bowtie $SAMPLE analyzed/.
	fi

	prnCmd "# FINISHED: COPYING TO REPO"

	# We run rsync again here to make sure everything copied and to
	# have the full transfer since we just changed $OUT with the
	# prnCmd above.
	if ! $DEBUG; then 
		rsync -avh --stats --exclude raw --exclude bowtie $SAMPLE analyzed/.
	fi
fi

#################################################
# Generate blast database from reads
#################################################
if [[ "$ACTION" == "blastdb" ]]; then
	prnCmd "# BEGIN: CREATE BLAST DATABASE"

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
	prnCmd "awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}' $SAMPLE/raw/unaligned_2.fq >> $SAMPLE/blast.db/raw.fa"
	if ! $DEBUG; then 
		awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $SAMPLE/raw/unaligned_2.fq >> $SAMPLE/blast.db/raw.fa
	fi

	# need to save current directory so we can return here. We also
	# need to adjust $OUT so prnCmd() still works when we change
	# directories
	prnCmd "CUR_DIR=`pwd`"
	CUR_DIR=`pwd`

	prnCmd "OUT_SAV=$OUT"
	OUT_SAV=$OUT

	prnCmd "cd $SAMPLE/blast.db"
	if ! $DEBUG; then 
		cd $SAMPLE/blast.db

		OUT=../../$OUT
		prnCmd "OUT=../../$OUT"
	fi

	prnCmd "makeblastdb -in raw.fa -dbtype nucl -out $SAMPLE -title \"$SAMPLE\""
	if ! $DEBUG; then 
		makeblastdb -in raw.fa -dbtype nucl -out $SAMPLE -title "$SAMPLE"
	fi

	prnCmd "rm raw.fa"
	rm raw.fa

	# return to proper directory and restore $OUT
	prnCmd "cd $CUR_DIR"
	if ! $DEBUG; then 
		cd $CUR_DIR

		OUT=$OUT_SAV
		prnCmd "OUT=$OUT_SAV"
	fi

	prnCmd "# FINISHED: CREATING BLAST DATABASE"
fi

#################################################
# Print out stats
#################################################
#if [[ "$ACTION" == "all" || "$ACTION" == "stats" ]]; then
if [[ "$ACTION" == "stats" ]]; then
	prnCmd "# BEGIN: STATS"

	echo "-- BLAST --"
	cat $SAMPLE/blast/species.txt; 

	NUMREADS=`grep "Total Lines" $SAMPLE/blast/species.txt | awk '{print $3}'`
	BLASTHITS=`grep "Number hits:" $SAMPLE/blast/species.txt | awk '{print $3}'`
	RAT=`grep "Rat Hits" $SAMPLE/blast/species.txt | awk '{print $3}'`
	MOUSE=`grep "Mouse Hits" $SAMPLE/blast/species.txt | awk '{print $3}'`
	HUMAN=`grep "Human Hits" $SAMPLE/blast/species.txt | awk '{print $3}'`
	BACTER=`grep "Bacteria Hits" $SAMPLE/blast/species.txt | awk '{print $3}'`
	
	#echo -e "\n-- FastQC --"
	#tar cvf $SAMPLE.fastqc.tar $SAMPLE/fastqc
	
	echo -e "\n-- Adapter Trimming --"
	cat $SAMPLE/trimAdapter.stats.txt; 
	
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
	echo -e "Num Reads\tNum Reads Trimmed\tRUM\tBlast Tot Hits\tRat\tMouse\tHuman\tBacteria\tBowtie\tAdapt Trim Kept\tAdapt Trim Disc\tPolyAT Trim Kept\tPolyAT Trim Disc\tDate\n"
	echo -e "$NUMREADS\t$NUMTRIM\t$RUMALI\t$BLASTHITS\t$RAT\t$MOUSE\t$HUMAN\t$BACTER\t\t$TRIMAD_K\t$TRIMAD_D\t$TRIMAT_K\t$TRIMAT_D\t$(date +%m/%d/%Y)\n"

	prnCmd "# FINISHED: STATS"
fi

exit 0

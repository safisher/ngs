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
# INPUT: $SAMPLE/init/unaligned_1.fq
# OUTPUT: $SAMPLE/blast/blast.txt (blast output), $SAMPLE/blast/species.txt (species hit counts)
# REQUIRES: blastn (provided with Blast version 2), randomSample.py, parseBlast.py
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` blast OPTIONS sampleID    --  run blast on randomly sampled subset of reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_BLAST() {
	echo -e "Usage:\n\t`basename $0` blast [-r numReads] [-k kmer] -p numProc -s species sampleID"
	echo -e "Input:\n\tsampleID/init/unaligned_1.fq"
	echo -e "Output:\n\tsampleID/blast/blast.txt (blast output)\n\tsampleID/blast/sampleID.blast.stats.txt (species hit counts)"
	echo -e "Requires:\n\tblastn ( ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ )\n\trandomSample.py ( https://github.com/safisher/ngs )\n\tparseBlast.py ( https://github.com/safisher/ngs )"
	echo -e "Options:"
	echo -e "\t-r numReads - number of reads to randomly select (default = 5000)"
	echo -e "\t-k kmer - check for presence of k-mer in reads that failed to align"
	echo -e "\t-p numProc - number of cpu to use"
	echo -e "\t-s species - expected species\n"
	echo -e "Run blast on 5000 reads randomly sampled from init/unaligned_1.fq. Blast paramters used are 'num_descriptions: 10 num_alignments: 10 word_size: 15 gapopen: 3 gapextend: 1 evalue: 1e-15'. The output is put in a directory called 'blast'. The species.txt file contains number of reads mapping to each species (mouse, rat, human, bacteria)."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

# number of reads to randomly select
ngsLocal_BLAST_NUM_READS=5000
ngsLocal_BLAST_KMER=""

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# BLAST args: -p value, sampleID
##########################################################################################

ngsArgs_BLAST() {
	if [ $# -lt 5 ]; then printHelp "BLAST"; fi
		
	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-p) NUMCPU=$2
				shift; shift;
				;;
			-s) SPECIES=$2
				shift; shift;
				;;
			-k) ngsLocal_BLAST_KMER=$2
				shift; shift;
				;;
			-r) ngsLocal_BLAST_NUM_READS=$2
				shift; shift;
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
# This will do a BLAST search on 5,000 untrimmed reads, using the nt database.
##########################################################################################

ngsCmd_BLAST() {
	prnCmd "# BEGIN: BLAST"

	# make relevant directory
	if [ ! -d $SAMPLE/blast ]; then 
		prnCmd "mkdir $SAMPLE/blast"
		if ! $DEBUG; then mkdir $SAMPLE/blast; fi
	fi
		
    # Get ngsLocal_BLAST_NUM_READS (5,000) randomly sampled reads
    # Usage: randomSample.py <num lines> <lines grouped> <input> <output>
	prnCmd "randomSample.py $ngsLocal_BLAST_NUM_READS 4 $SAMPLE/init/unaligned_1.fq $SAMPLE/blast/raw.fq > $SAMPLE/blast/sampling.out.txt"
	if ! $DEBUG; then 
		randomSample.py $ngsLocal_BLAST_NUM_READS 4 $SAMPLE/init/unaligned_1.fq $SAMPLE/blast/raw.fq > $SAMPLE/blast/sampling.out.txt
	fi
	
    # Convert fastq file to fasta file
	prnCmd "awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}' $SAMPLE/blast/raw.fq > $SAMPLE/blast/raw.fa"
	if ! $DEBUG; then 
		awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $SAMPLE/blast/raw.fq > $SAMPLE/blast/raw.fa
	fi
	
    # Run BLAST. Output file should end with ".txt"
	prnCmd "blastn -query $SAMPLE/blast/raw.fa -db nt -num_descriptions 10 -num_alignments 10 -word_size 15 -gapopen 3 -gapextend 1 -evalue 1e-15 -num_threads $NUMCPU -out $SAMPLE/blast/blast.txt"
	if ! $DEBUG; then 
		blastn -query $SAMPLE/blast/raw.fa -db nt -num_descriptions 10 -num_alignments 10 -word_size 15 -gapopen 3 -gapextend 1 -evalue 1e-15 -num_threads $NUMCPU -out $SAMPLE/blast/blast.txt
	fi
	
    # Parse BLAST output. Will generate *.cvs and *.hits files.
    # Usage: parseBlast.py targetSpecies readsFastaFile blastFile
	prnCmd "parseBlast.py $SPECIES $SAMPLE/blast/raw.fa $SAMPLE/blast/blast.txt"
	if ! $DEBUG; then 
		parseBlast.py $SPECIES $SAMPLE/blast/raw.fa $SAMPLE/blast/blast.txt

		# test for kmer in reads that didn't map with BLAST. These reads
		# are put into the file "failed.fa" by parseBlast.py. 
		if [[ ! -z "$ngsLocal_BLAST_KMER" ]]; then
			# $ngsLocal_BLAST_KMER has a non-zero length, so user gave us
			# a k-mer to test
			local count=$(grep -c "$ngsLocal_BLAST_KMER" $SAMPLE/blast/failed.fa)
			
			# We add the count to THE TOP of the speciesCounts.txt file
			# that parseBlast.py generates. We use the top of the file
			# because the last line in the file contains a tab-delimited
			# string of counts generated by parseBlast.py.
			local line="k-mer count: ${count}"
			if [[ ${OS_VERSION} == "Darwin" ]]; then
				# Sed on the Mac requires a zero-length extension when
				# using the -i option. This is not a requirement in linux.
				sed -i "" '1 i\
'"$line"'
' $SAMPLE/blast/speciesCounts.txt
			else
				sed -i "1 i${line}" $SAMPLE/blast/speciesCounts.txt
			fi
			
			line="k-mer sequence: ${ngsLocal_BLAST_KMER}"
			if [[ ${OS_VERSION} == "Darwin" ]]; then
				sed -i "" '1 i\
'"$line"'
' $SAMPLE/blast/speciesCounts.txt
			else
				sed -i "1 i${line}" $SAMPLE/blast/speciesCounts.txt
			fi
		fi
	fi

	# rename output stats file to conform to other modules
	prnCmd "mv $SAMPLE/blast/speciesCounts.txt $SAMPLE/blast/$SAMPLE.blast.stats.txt"
	if ! $DEBUG; then 
		mv $SAMPLE/blast/speciesCounts.txt $SAMPLE/blast/$SAMPLE.blast.stats.txt
	fi
	
	# run error checking
	if ! $DEBUG; then ngsErrorChk_BLAST $@; fi

    # print version info in $SAMPLE directory. We do this AFTER
    # parseBlast.py has run because we need to get the version number
    # from the output file.
	prnCmd "# blastn version: blastn -version | tail -1 | awk '{print \$3}' | sed s/,//"
	prnCmd "# parseBlast.py version: grep parseBlast $SAMPLE/blast/$SAMPLE.blast.stats.txt | awk -F: '{print \$2}'"
	if ! $DEBUG; then 
		# gets this: "Package: blast 2.2.28, build Mar 12 2013 16:52:31"
		# returns this: "2.2.28"
		ver=$(blastn -version | tail -1 | awk '{print $3}' | sed s/,//)
		ver1=$(grep parseBlast $SAMPLE/blast/$SAMPLE.blast.stats.txt | awk -F: '{print $2}')
		prnVersion "blast" "program\tversion\tprogram\tversion\tspecies" "blastn\t$ver\tparseBlast.py\t$ver1\t$SPECIES"
	fi

	prnCmd "# FINISHED: BLAST"
}

##########################################################################################
# ERROR CHECKING. Make sure output file exists and contains species counts.
##########################################################################################

ngsErrorChk_BLAST() {
	prnCmd "# BLAST ERROR CHECKING: RUNNING"

	inputFile="$SAMPLE/init/unaligned_1.fq"
	outputFile="$SAMPLE/blast/$SAMPLE.blast.stats.txt"

	# make sure expected output file exists
	if [ ! -f $outputFile ]; then
		errorMsg="Expected BLAST output file does not exist.\n"
		errorMsg+="\tinput file: $inputFile\n"
		errorMsg+="\toutput file: $outputFile\n"
		prnError "$errorMsg"
	fi

	# compute number of lines in speciesCounts.txt
	counts=`wc -l $outputFile | awk '{print $1}'`

	# if counts file has less than 3 lines, then BLAST didn't work
	if [ "$counts" -lt "3" ]; then
		errorMsg="BLAST failed to run properly and there are no species counts.\n"
		errorMsg+="\tinput file: $inputFile\n"
		errorMsg+="\toutput file: $outputFile\n"
		prnError "$errorMsg"
	fi

	prnCmd "# BLAST ERROR CHECKING: DONE"
}

##########################################################################################
# PRINT STATS. Prints a tab-delimited list stats of interest.
##########################################################################################

ngsStats_BLAST() {
	if [ $# -ne 1 ]; then
		prnError "Incorrect number of parameters for ngsStats_BLAST()."
	fi
	
	local statsFile="$SAMPLE/blast/$SAMPLE.blast.stats.txt"

	# the second to the last line of the stats file is a tab-delimited lists of headers
	# Total Hits	Hits Not Counted	Bacteria	Fish	Fly	Human	Mouse	Rat	Yeast
	local header=$(tail -2 $statsFile | head -1)
	# the last line of the stats file is a tab-delimited lists of values
	local values=$(tail -1 $statsFile)

	# test if k-mer search results exist in stats file. If so the
	# include in output.
	local kmerSeq=$(grep "k-mer sequence" $statsFile | awk '{print $3}')
	# if not empty, then k-mer results exist
	if [[ ! -z "$kmerSeq" ]]; then
		# append sequence to header/values list
		header+="\t$kmerSeq"

		# add the count to the values list
		local kmerCount=$(grep "k-mer count" $statsFile | awk '{print $3}')
		# append sequence to header/values list
		values+="\t$kmerCount"
	fi

	case $1 in
		header)
			echo $header
			;;

		values) 
			echo $values
			;;

		keyvalue) 
			# output key:value pair of stats

			# the bash IFS variable dictates the word delimiting which is " \t\n" 
			# by default. We want to only delimite by tabs for the case here.
			local IFS=$'\t'

			# convert tab-delimited header/values variables to array
			declare -a headerArray=($header)
			declare -a valuesArray=($values)
			
			# output a tab-delimited, key:value list
			numFields=${#headerArray[@]}
			for ((i=0; i<$numFields-1; ++i)); do
				echo -en "${headerArray[$i]}:${valuesArray[$i]}\t"
			done
			echo "${headerArray[$numFields-1]}:${valuesArray[$numFields-1]}"
			;;

		*) 
			# incorrect argument
			prnError "Invalid parameter for ngsStats_BLAST() (got $1, expected: 'header|values')."
			;;
	esac
}

#!/bin/bash

# Copyright (c) 2012,2013,2015 Stephen Fisher, Jamie Shallcross, and Junhyong Kim, 
# University of Pennsylvania.  All Rights Reserved.
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
# INPUT: $SAMPLE/star/STAR_Unique.bam
# OUTPUT: $SAMPLE/verse/$SAMPLE.verse.cnts.txt, $SAMPLE/verse/$SAMPLE.verse.log.txt, $SAMPLE/verse/$SAMPLE.verse.err.txt
# REQUIRES: VERSE
##########################################################################################

## Version 0.2.0

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` verse OPTIONS sampleID    --  run VERSE on unique mapped reads\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_VERSE() {
    echo -e "Usage:\n\t`basename $0` verse [-i inputDir] [-f inputFile] [-stranded] [-l features] [-lines_sines] [-id idAttr] -p numProc -s species sampleID"
    echo -e "Input:\n\tsampleID/inputDir/inputFile"
    echo -e "Output:\n\tsampleID/verse/sampleID.verse.cnts.txt\n\tsampleID/verse/sampleID.verse.log.txt\n\tsampleID/verse/sampleID.verse.err.txt"
    echo -e "Requires:\n\tVERSE version 0.1.1 or later\n"
    echo -e "Options:"
    echo -e "\t-i inputDir - location of source file (default: star)."
    echo -e "\t-f inputFile - source file (default: sampleID.star.unique.bam)."
    echo -e "\t-stranded - use strand information (default: no)."
    echo -e "\t-p numProc - maximum number of cpu to use."
    echo -e "\t-l features - ordered, comma or semicolon separated list of feature types for verse hierarchical assignment. (default: exon)"
    echo -e "\t-lines_sines - also compute line and sine counts, using intersection-nonempty (default: no) This count will be independant of counts for exons/introns/intergenic regions. The output file will have the suffix 'xine'."
    echo -e "\t-id idAttr - name of the the GTF field that contains the name/ID of a gene (default: gene_id). Counts will sum all features with the same value for this field."
    echo -e "\t-s species - species from repository: $verse_REPO.\n"
    echo -e "Run VERSE gene quantification. This requires a BAM file as generated by either RUMALIGN or STAR (STAR by default)."
    echo -e "The following VERSE parameter values are used for exon counting:\n \t-z 3 (intersection_nonempty) "
    echo -e "Each type of feature being counted (i.e. exons, introns, intergentic, mitochondrial, lines and sines) will be run in sucession, and generate separate counts files, in the form:"
    echo -e "\tSampleID.verse.exon.cnts.txt: exon counts"
    echo -e "\tSampleID.verse.intron.cnts.txt: intron counts"
    echo -e "For a description of the VERSE parameters see [[VERSE documentation]]\n"
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_verse_INP_DIR="star"
# the default for ngsLocal_verse_INP_FILE is set in ngsCmd_VERSE()
# because it depends on the value of $SAMPLE and $SAMPLE doesn't have
# a value until the ngsCmd_VERSE() function is run.
ngsLocal_verse_INP_FILE=""
ngsLocal_verse_STRANDED="0"
ngsLocal_verse_EXONS="1"
ngsLocal_verse_MITO="0"
ngsLocal_verse_INTRONS="0"
ngsLocal_verse_INTERGENIC="0"
ngsLocal_verse_xINEs="0"
ngsLocal_verse_FEATURES="exon"
ngsLocal_VERSE_NUMCPU=1


# use "gene_id" by default but let users change to "gene_name" or whatever needed
ngsLocal_verse_ID_ATTR="gene_id"

# Run verse in intersection-nonempty mode 
ngsLocal_verse_MODE="3"

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# verse args: -s value, -g value, sampleID
##########################################################################################

ngsArgs_VERSE() {
    if [ $# -lt 3 ]; then printHelp "VERSE"; fi
    
    # getopts doesn't allow for optional arguments so handle them manually
    while true; do
	if [[ -z $1 ]] ; then printHelp "VERSE"; fi
	case $1 in
	    -i) ngsLocal_verse_INP_DIR=$2
		shift; shift;
		;;
	    -f) ngsLocal_verse_INP_FILE=$2
		shift; shift;
		;;
	    -p) ngsLocal_VERSE_NUMCPU=$2
		shift; shift;
		;;
	    -l) ngsLocal_verse_FEATURES=$2
		shift; shift;
		;;
	    -stranded) ngsLocal_verse_STRANDED="1"
		shift;
		;;
	    -lines_sines) ngsLocal_verse_xINEs="1"
		shift;
		;;
	    -id) ngsLocal_verse_ID_ATTR=$2
		shift; shift;
		;;
	    -s) SPECIES=$2
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
# Run VERSE on uniqely mapped alignments, as generated by the POST command.
##########################################################################################

ngsCmd_VERSE() { 
    prnCmd "# BEGIN: verse"
    
    if [[ $ngsLocal_verse_FEATURES == *"exon"* ]]; then
	ngsLocal_verse_EXONS="1"
    fi

    if [[ $ngsLocal_verse_FEATURES == *"mito"* ]]; then
	ngsLocal_verse_MITO="1" 
    fi

    if [[ $ngsLocal_verse_FEATURES == *"intron"* ]]; then
	ngsLocal_verse_INTRONS="1" 
    fi

    if [[ $ngsLocal_verse_FEATURES == *"intergenic"* ]]; then
	ngsLocal_verse_INTERGENIC="1"
    fi 

    # make relevant directory
    if [ ! -d $SAMPLE/verse ]; then 
	prnCmd "mkdir $SAMPLE/verse"
	if ! $DEBUG; then mkdir $SAMPLE/verse; fi
    fi
    
    # print version info in $SAMPLE directory
    prnCmd "# VERSE version: verse -v"
    if ! $DEBUG; then 
	# returns: "VERSE v0.1.1" or similar
	ver=$(verse -v 2>&1 | grep 'VERSE' | cut -d' ' -f2)
	prnVersion "verse" \
	"program\tversion\ttranscriptome\tstranded\tID_attribute\tintrons\tintergenic\tlines-sines" \
	"verse\t$ver\t${verse_REPO}/${SPECIES}.gtf\t$ngsLocal_verse_STRANDED\t$ngsLocal_verse_ID_ATTR\t$ngsLocal_verse_INTRONS\t$ngsLocal_verse_INTERGENIC\t$ngsLocal_verse_xINEs"
    fi
    
    # if the user didn't provide an input file then set it to the
    # default
    if [[ -z $ngsLocal_verse_INP_FILE ]]; then 
	ngsLocal_verse_INP_FILE="$SAMPLE/star/$SAMPLE.star.unique.bam"
    fi 
    # We assume that the alignment file exists



    ngsLocal_verse_FEATURES=`sed s/,/\;/g <<< $ngsLocal_verse_FEATURES`

    local THREADS=$(( $ngsLocal_VERSE_NUMCPU > 5 ? 2 : 1)) # Testing shows no benefit from >2 helper threads in verse

    
    prnCmd "verse -T $THREADS -t \"$ngsLocal_verse_FEATURES\" -z $ngsLocal_verse_MODE --nonemptyModified -g $ngsLocal_verse_ID_ATTR -s $ngsLocal_verse_STRANDED -R -a $verse_REPO/$SPECIES.gtf -o $SAMPLE/verse/$SAMPLE.verse $ngsLocal_verse_INP_FILE > $SAMPLE/verse/$SAMPLE.verse.out" 

    verse -T $THREADS -t "$ngsLocal_verse_FEATURES" -z $ngsLocal_verse_MODE --nonemptyModified -g $ngsLocal_verse_ID_ATTR -s $ngsLocal_verse_STRANDED -R -a "$verse_REPO/$SPECIES.gtf" -o $SAMPLE/verse/$SAMPLE.verse $ngsLocal_verse_INP_FILE > $SAMPLE/verse/$SAMPLE.verse.out &
    
    if [[ $ngsLocal_VERSE_NUMCPU -le 3 ]]; then
	wait
    fi

    if [ $ngsLocal_verse_xINEs = "1" ]; then
	prnCmd "verse -T $THREADS -t "xine" -z $ngsLocal_verse_MODE --nonemptyModified -g $ngsLocal_verse_ID_ATTR -s $ngsLocal_verse_STRANDED -a "$verse_REPO/$SPECIES.gtf" -o $SAMPLE/verse/$SAMPLE.verse.lines_sines $ngsLocal_verse_INP_FILE > $SAMPLE/verse/$SAMPLE.verse.lines_sines.out &"
	
	verse -T $THREADS -t "xine" -z $ngsLocal_verse_MODE --nonemptyModified -g $ngsLocal_verse_ID_ATTR -s $ngsLocal_verse_STRANDED -a "$verse_REPO/$SPECIES.gtf" -o $SAMPLE/verse/$SAMPLE.verse.lines_sines $ngsLocal_verse_INP_FILE > $SAMPLE/verse/$SAMPLE.verse.lines_sines.log.txt &
    fi

    wait


    prnCmd "# splitting log files and sorting output counts"
    # Postprocess output files. Sort gene counts, produce warnings/errors, and gzip detail file. 
    if ! $DEBUG; then 
	# only generate error file if Warnings exist. If we run grep
	# and it doesn't find any matches then it will exit with an
	# error code which would cause the program to crash since we
	# use "set -o errexit"
	ngsPostverse
    fi

    # run error checking
    if ! $DEBUG; then 

	#Using a subshell since we're changing IFS temporarily 
	( # BEGIN SUBSHELL
	    IFS=';'
	    for feature in $ngsLocal_verse_FEATURES; do
		ngsErrorChk_VERSE $feature
	    done
	) # END SUBSHELL

	if [ $ngsLocal_verse_xINEs = "1" ]; then ngsErrorChk_VERSE lines_sines; fi

    fi
    
    prnCmd "# FINISHED: verse"
}

ngsPostverse() {

    # only generate error file if Warnings exist. If we run grep
    # and it doesn't find any matches then it will exit with an
    # error code which would cause the program to crash since we
    # use "set -o errexit"
    local containsWarningsI=$(grep -ci 'WARNING' $SAMPLE/verse/$SAMPLE.verse.out)
    if [[ $containsWarningsI -gt 0 ]]; then
	prnCmd "grep -i 'WARNING' $SAMPLE/verse/$SAMPLE.verse.out > $SAMPLE/verse/$SAMPLE.verse.err.txt"
	grep -i 'WARNING' $SAMPLE/verse/$SAMPLE.verse.out > $SAMPLE/verse/$SAMPLE.verse.err.txt
    fi

    prnCmd "mv $SAMPLE/verse/$SAMPLE.verse.out $SAMPLE/verse/$SAMPLE.verse.log.txt"
    mv $SAMPLE/verse/$SAMPLE.verse.out $SAMPLE/verse/$SAMPLE.verse.log.txt

    if [[ -f $SAMPLE/verse/$SAMPLE.verse.exon.summary.txt ]]; then # If only exons are counted then verse writes to a different name, which messes up stats. 
	prnCmd "mv $SAMPLE/verse/$SAMPLE.verse.exon.summary.txt $SAMPLE/verse/$SAMPLE.verse.summary.txt"
	mv $SAMPLE/verse/$SAMPLE.verse.exon.summary.txt $SAMPLE/verse/$SAMPLE.verse.summary.txt
    fi
    
    #Using a subshell since we're changing IFS temporarily 
    ( # BEGIN SUBSHELL
	IFS=';'
	for feature in $ngsLocal_verse_FEATURES; do
	    #The subshell is used so that we can sort everything but the first line which is a header
	    prnCmd "(head -1 $SAMPLE/verse/$SAMPLE.verse.${feature}.txt && tail -n +2 $SAMPLE/verse/$SAMPLE.verse.${feature}.txt | sort) > $SAMPLE/verse/$SAMPLE.verse.${feature}.cnts.txt"
	    (head -1 $SAMPLE/verse/$SAMPLE.verse.${feature}.txt && tail -n +2 $SAMPLE/verse/$SAMPLE.verse.${feature}.txt | sort -V) > $SAMPLE/verse/$SAMPLE.verse.${feature}.cnts.txt

	    prnCmd "rm $SAMPLE/verse/$SAMPLE.verse.${feature}.txt"
	    rm $SAMPLE/verse/$SAMPLE.verse.${feature}.txt
	done
    ) # END SUBSHELL


    if [[ -f $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.txt ]]; then
	prnCmd "(head -1 $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.txt && tail -n +2 $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.txt | sort) > $SAMPLE/verse/$SAMPLE.verse.lines_sines.cnts.txt"
	(head -1 $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.txt && tail -n +2 $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.txt | sort -V) > $SAMPLE/verse/$SAMPLE.verse.lines_sines.cnts.txt

	prnCmd "rm $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.txt"
	rm $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.txt

	prnCmd "mv $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.summary.txt $SAMPLE/verse/$SAMPLE.verse.lines_sines.summary.txt"
	mv $SAMPLE/verse/$SAMPLE.verse.lines_sines.xine.summary.txt $SAMPLE/verse/$SAMPLE.verse.lines_sines.summary.txt
    fi
    
    prnCmd "gzip -f $SAMPLE/verse/$SAMPLE.verse.detail.txt"
    gzip -f $SAMPLE/verse/$SAMPLE.verse.detail.txt    
}


##########################################################################################
# ERROR CHECKING. Make sure output file exists, is not effectively
# empty and warn user if VERSE output any warnings.
##########################################################################################

# $1 should be "exons" or "introns" or "mito"
ngsErrorChk_VERSE() {
    prnCmd "# verse ERROR CHECKING $1: RUNNING"

    inputFile="$SAMPLE/$ngsLocal_verse_INP_DIR/$ngsLocal_verse_INP_FILE"
    outputFile="$SAMPLE/verse/$SAMPLE.verse.${1}.cnts.txt"
    
    # make sure expected output file exists
    if [ ! -f $outputFile ]; then
	errorMsg="Expected VERSE output file does not exist.\n"
	errorMsg+="\tinput file: $inputFile\n"
	errorMsg+="\toutput file: $outputFile\n"
	prnError "$errorMsg"
    fi

    # if cnts file only has 1 line then error and print contents of log file
    counts=`wc -l $outputFile | awk '{print $1}'`
    # if counts file only has one line, then VERSE didn't work
    if [ "$counts" -eq "1" ]; then
	errormsg="verse failed to run properly. see verse error below:\n"
	errormsg+="\tinput file: $inputFile\n"
	errormsg+="\toutput file: $outputFile\n\n"
	errormsg+=`cat $SAMPLE/verse/$SAMPLE.verse.log.txt`
	prnerror "$errormsg"
    fi
    
    # Check err file for errors
    if [ -s $SAMPLE/verse/$SAMPLE.verse.${1}err.txt ]; then
	warningMsg="Review the error file listed below to view VERSE warnings.\n"
	warningMsg+="\tinput file: $inputFile\n"
	warningMsg+="\toutput file: $outputFile\n"
	warningMsg+="\tERROR FILE: $SAMPLE/verse/$SAMPLE.verse.${1}.err.txt\n"
	prnWarning "$warningMsg"
    fi
    
    prnCmd "# verse ERROR CHECKING $1: DONE"
}

##########################################################################################
# PRINT STATS. Prints a tab-delimited list stats of interest.
##########################################################################################

ngsStats_VERSE() {
    if [ $# -ne 1 ]; then
	prnError "Incorrect number of parameters for ngsStats_VERSE()."
    fi
    

    noFEATURE=$(grep -iw "NoFeature" $SAMPLE/verse/$SAMPLE.verse.summary.txt | awk '{print $2}')
    missingMATES=$(grep -iw "MissingMates" $SAMPLE/verse/$SAMPLE.verse.summary.txt | awk '{print $2}')
    totalReads=$(grep -iw "TotalReadPairs" $SAMPLE/verse/$SAMPLE.verse.summary.txt | awk '{print $2}')
    
    header=""
    values=""
    percsHeader=""
    percsValues=""

    if [ -f $SAMPLE/verse/${SAMPLE}.verse.exon.cnts.txt ]; then 
	ngsHelperStatsverse "exon" "exons Level 1,2"
	#else
	#    header="$header\t\t\t\t\t\t"
	#    values="$values\t\t\t\t\t\t"

	pExonL1L2=$(grep AssignedExonFraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%')
	#spikeReadsCounted=$(grep "spikeIn" $SAMPLE/verse/${SAMPLE}.verse.exon.cnts.txt | awk -F '\t' '{sum += $2} END {print sum}')
	spikeReadsCounted=$(grep "spikeIn" $SAMPLE/verse/${SAMPLE}.verse.exon.cnts.txt | awk -F '\t' '{sum += $2} END {print sum}')
	if [[ -z $spikeReadsCounted ]]; then
	    spikePerc=0.0
	else
	    spikePerc=$(echo $spikeReadsCounted $totalReads | awk -F' ' '{printf "%.2f", $1/$2 * 100}' )
	fi
	pExonL1L2=$(echo $pExonL1L2 - $spikePerc | bc)
	percsHeader="$percsHeader\tPerc: Exons Level 1,2"
	percsValues="$percsValues\t$pExonL1L2"
    fi

    if [ -f $SAMPLE/verse/${SAMPLE}.verse.exon-lev3.cnts.txt ]; then 
	ngsHelperStatsverse "exon-lev3" "exons Level 3"
	#else
	#    header="$header\t\t\t\t\t\t"
	#    values="$values\t\t\t\t\t\t"
	pExonL3=`grep AssignedExon-lev3Fraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Exons Level 3"
	percsValues="$percsValues\t$pExonL3"
    fi

    if [ -f $SAMPLE/verse/${SAMPLE}.verse.anti-exon.cnts.txt ]; then 
	ngsHelperStatsverse "anti-exon" "anti-exons"
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
	pAntiExon=`grep AssignedAnti-exonFraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Anti-exons"
	percsValues="$percsValues\t$pAntiExon"
    fi

    if [ -f $SAMPLE/verse/${SAMPLE}.verse.intron-lev1-lev2.cnts.txt ]; then 
	ngsHelperStatsverse "intron-lev1-lev2" "Introns Level 1,2"
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
	pIntron=`grep AssignedIntron-lev1-lev2Fraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Introns Level 1,2"
	percsValues="$percsValues\t$pIntron"
    fi

    if [ -f $SAMPLE/verse/${SAMPLE}.verse.intron-lev3.cnts.txt ]; then 
	ngsHelperStatsverse "intron-lev3" "Introns Level 3"
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
	pL3Intron=`grep AssignedIntron-lev3Fraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Introns Level 3"
	percsValues="$percsValues\t$pL3Intron"
    fi

    if [ -f $SAMPLE/verse/${SAMPLE}.verse.intron.cnts.txt ]; then 
	ngsHelperStatsverse "intron" "Introns"
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
	pIntron=`grep AssignedIntronFraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Introns"
	percsValues="$percsValues\t$pIntron"
    fi

    if [ -f $SAMPLE/verse/${SAMPLE}.verse.anti-intron.cnts.txt ]; then 
	ngsHelperStatsverse "anti-intron" "Anti-introns"
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
	pAntiIntron=`grep AssignedAnti-intronFraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Anti-introns"
	percsValues="$percsValues\t$pAntiIntron"
    fi

    if [ -f $SAMPLE/verse/$SAMPLE.verse.mito.cnts.txt ]; then 
	ngsHelperStatsverse "mito" "Mito"
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
	pMito=`grep AssignedMitoFraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Mito"
	percsValues="$percsValues\t$pMito"
    fi

    if [ -f $SAMPLE/verse/$SAMPLE.verse.anti-mito.cnts.txt ]; then 
	ngsHelperStatsverse "anti-mito" "Anti-mito"
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
	pAntiMito=`grep AssignedAnti-mitoFraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Anti-mito"
	percsValues="$percsValues\t$pAntiMito"
    fi

    if [ -f $SAMPLE/verse/$SAMPLE.verse.intergenic.cnts.txt ]; then 
	ngsHelperStatsverse "intergenic" "Intergenic"
	#else
	#    header="$header\t\t\t\t"
	#    values="$values\t\t\t\t"
	pIntergenic=`grep AssignedIntergenicFraction $SAMPLE/verse/${SAMPLE}.verse.summary.txt | cut -f2 | tr -d '%'`
	percsHeader="$percsHeader\tPerc: Intergenic"
	percsValues="$percsValues\t$pIntergenic"
    fi

    if [ -f $SAMPLE/verse/$SAMPLE.verse.lines_sines.cnts.txt ]; then 
	ngsHelperStatsverse_LS
	#else
	#    header="$header\t\t\t\t\t"
	#    values="$values\t\t\t\t\t"
    fi

    header="$header\tNo Feature\tReads Missing Mates${percsHeader}\tPerc: Spike-In"
    values="$values\t$noFEATURE\t$missingMATES$percsValues\t$spikePerc"

    case $1 in
	header) 
	    echo "$header"
	    ;;
	
	values) 
	    echo "$values"
	    ;;
	
	*) 
	    # incorrect argument
	    prnError "Invalid parameter for ngsStats_VERSE() (got $1, expected: 'header|values')."
	    ;;
    esac
}

ngsHelperStatsverse() {
    # $1 = exons | introns | mito | intergenic
    # $2 = ie "gene (exons) or "intergenic region"

    local gene="Gene"
    if [[ $1 == "intergenic" ]]; then gene="Region"; fi

    if [[ $1 == "exon" ]]; then
	# total number of reads that mapped unambigously to genes
	readsCounted=$(grep -v "spikeIn." $SAMPLE/verse/$SAMPLE.verse.${1}.cnts.txt | awk -F '\t' '{sum += $2} END {print sum}')
	header="${header}$2: non-spikeIn Reads Counted"
	values="${values}$readsCounted"
	
	# total number of reads that mapped unambigously to ERCC controls
	spikeReadsCounted=$(grep "spikeIn." $SAMPLE/verse/$SAMPLE.verse.${1}.cnts.txt | awk -F '\t' '{sum += $2} END {print sum}')
	header="$header\t$2: spikeIn Reads Counted"
	values="$values\t$spikeReadsCounted"
    else
	readsCounted=$($GREPP -i "^Assigned${1}\s" $SAMPLE/verse/$SAMPLE.verse.summary.txt | awk '{print $2}')
	header="$header\t$2: Reads Counted"
	values="$values\t$readsCounted"
    fi

    # number of genes with at least 1 read mapped
    numGenes=$($GREPP -v "\t0$" $SAMPLE/verse/$SAMPLE.verse.${1}.cnts.txt | grep -v "gene" | wc -l)
    header="$header\t$2: Num ${gene}s"
    values="$values\t$numGenes"

    
    # average number of reads that mapped unambigously to genes
    if [[ $numGenes -gt 0 ]]; then
	avgReadPerGene=$(($readsCounted/$numGenes))
    else
	avgReadPerGene=0
    fi
    header="$header\t$2: Avg Read Per $gene"
    values="$values\t$avgReadPerGene"
    
    # maximum number of reads that mapped unambigously to a single gene
    maxReadsPerGene=$(grep -v "gene" $SAMPLE/verse/$SAMPLE.verse.${1}.cnts.txt | awk -F '\t' '{if(max=="") {max=$2}; if($2>max) {max=$2};} END {print max}')
    header="$header\t$2: Max Reads Per $gene"
    values="$values\t$maxReadsPerGene"
    
    # number of reads that didn't map to a gene region
    #noFeature=$(tail -9 $SAMPLE/verse/$SAMPLE.verse.summary.txt | head -1 | awk '{print $2}')
    #header="$header\tNo Feature (${1})"
    #values="$values\t$noFeature"
    
    # number of reads that completely overlapped two or more gene regions
    ambiguousMapped=$($GREPP -i "^Ambiguous${1}\s" $SAMPLE/verse/$SAMPLE.verse.summary.txt | awk '{print $2}')
    header="$header\t$2: Ambiguous Mapped"
    values="$values\t$ambiguousMapped"
}

ngsHelperStatsverse_LS() {
	# different function for lines & sines	

	# total number of reads that mapped unambigously to genes
	LINEreadsCounted=$(grep 'LINE'  $SAMPLE/verse/$SAMPLE.verse.lines_sines.cnts.txt | awk -F '\t' '{sum += $2} END {print sum}')
	header="$header\tLINE Reads"
	values="$values\t$LINEreadsCounted"
	
	SINEreadsCounted=$(grep 'SINE'  $SAMPLE/verse/$SAMPLE.verse.lines_sines.cnts.txt | awk -F '\t' '{sum += $2} END {print sum}')
	header="$header\tSINE Reads"
	values="$values\t$SINEreadsCounted"

	# number of genes with at least 1 read mapped
	numLINEs=$($GREPP -v "\t0$" $SAMPLE/verse/$SAMPLE.verse.lines_sines.cnts.txt | grep -v "gene" | grep 'LINE' | wc -l)
	header="$header\tNum LINEs"
	values="$values\t$numLINEs"
	
	# number of genes with at least 1 read mapped
	numSINEs=$($GREPP -v "\t0$" $SAMPLE/verse/$SAMPLE.verse.lines_sines.cnts.txt | grep -v "gene" | grep 'SINE' | wc -l)
	header="$header\tNum SINEs"
	values="$values\t$numSINEs"

	# number of reads that didn't map to a gene region
	noFeature=$(grep "NoFeature" $SAMPLE/verse/$SAMPLE.verse.lines_sines.summary.txt | awk '{print $2}')
	header="$header\tNeither LINEs nor SINEs"
	values="$values\t$noFeature"
}

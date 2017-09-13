#!/bin/bash

# Copyright (c) 2015, Jamie Shallcross and Junhyong Kim, University of
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

set -u
VERSION=0.1
SYSTEM="kimclust"

#from https://coderwall.com/p/q-ovnw/killing-all-child-processes-in-a-shell-script
#kills pipeline instances started in the background which would otherwise be messily orphaned
die() {
    # Get process group id
    PGID=$(ps -o pgid= $$ | grep -o [0-9]*)

    # Kill it in a new new process group
    setsid kill -- -$PGID
    exit 0
}

trap "die" SIGINT SIGTERM


if [[ $# -lt 1 || $1 == "-h" ]]; then 
    echo "launcher.sh [config_file] [threads]"
    echo
    echo "Launcher pipelines with appropriate settings for a group of samples."
    echo "Pipelines will either be bsub'd or started in the background, depending on system." 
    echo
    echo "config_file should be in the format-"
    printf "\tSampleID,SeqType,Species,SE/PE,Length,Stranded?,SpikeIn,Contaminants,TISA?\n"
    exit
fi

# Example set of configurations for parameters and features of interest for a given species
declare -A genomes=(["human"]="-s hg38.gencode25 -id gene_name -f exon,anti-exon,intron,anti-intron,intergenic -lines_sines")
genomes+=(["mouse"]="-s mm10 -id gene_id -f exon,intron,intergenic -lines_sines")
genomes+=(["rat"]="-s rn6 -id gene_id -f exon,intron,intergenic -lines_sines") 
genomes+=(["zebrafish"]="-s GRCz10 -id gene_name -f exon,intron")
genomes+=(["C. elegans"]="-s cel235 -id gene_name -f exon")

# Translate protocol name to contaminants file name
declare -A contaminants=(["Default"]="contaminants.fa")
contaminants+=(["Nextera"]="nextera.fa")
contaminants+=(["Smartseq"]="smartseq.fa")

# Default number of cpu is 6
if [[ -n $2 ]]; then 
    cpu=$2
else
    cpu=6
fi

while read line
do
    if [[ $line == \#* ]] #skip commented out lines
	then continue 
    fi

    FLAGS=""
    echo
    echo $line
    if [[ $line == "" ]]; then break; fi
    SAMPLE="Sample_`cut -d',' -f1 <<< $line`"
    SEQ_TYPE=`cut -d',' -f2 <<< $line`
    SPECIES=`cut -d',' -f3 <<< $line`
    ENDED=`cut -d',' -f4 <<< $line`
    LENGTH=`cut -d',' -f5 <<< $line`
    STRANDED=`cut -d',' -f6 <<< $line`
    SPIKE=`cut -d',' -f7 <<< $line` # Not currently used
    CONTAM=`cut -d',' -f8 <<< $line` 
    BARCODE=`cut -d',' -f9 <<< $line`
    NONSEQ_BC=`cut -d',' -f10 <<< $line`

    echo $CONTAM
    
    if [[ $ENDED == "SE" ]]; then FLAGS="$FLAGS -se"; fi
    if [[ $STRANDED == "TRUE" ]]; then FLAGS="$FLAGS -stranded"; fi


    FLAGS="$FLAGS -l $LENGTH"
    FLAGS="$FLAGS -c ${contaminants[$CONTAM]}"
    FLAGS="$FLAGS ${genomes[$SPECIES]}"

    if [[ $BARCODE == "TRUE" || $BARCODE == "True" ]]; then 
	SEQ_TYPE="RNASeq_TISA"
	FLAGS="$FLAGS -b $NONSEQ_BC"
    elif [[ $BARCODE == "APRA" ]]; then
	( apra; ) &
	#FLAGS="$FLAGS -noinit"
	continue
    fi

    echo "ngs.sh pipeline -t $SEQ_TYPE $FLAGS -chgrp repo-admin -p $cpu $SAMPLE &"
    ngs.sh pipeline -t $SEQ_TYPE $FLAGS -chgrp repo-admin -p $cpu $SAMPLE &
    sleep 5 #Stop it from bunching up all the IO

#done
done < $1

wait

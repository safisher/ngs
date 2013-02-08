#!/bin/bash

# Copyright (c) 2011,2012,2013, Stephen Fisher and Junhyong Kim, University of
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


if [ -z "$4" ]; then
	echo "Paired alignment with bowtie, followed by aligning the unmapped 1st and 2nd sequences."
	echo "Unaligned data is expected to be named 'unaligned_1.fq' and 'unaligned_2.fq'"
	echo "runBowtiePaired.sh <input dir> <output dir> <ref seq -> hg19 | mm9 | rat.m4> <num nodes>"
	exit 0 
fi

INP_DIR=$1
OUT_DIR=$2
REF_SEQ="/lab/repo/resources/indexes/$3_genome"
NODES=$4

OUT_FILE="$OUT_DIR/stats.txt"

READ_1=$INP_DIR/unaligned_1.fq
READ_2=$INP_DIR/unaligned_2.fq

echo -e "Starting bowtie..." > $OUT_FILE
date >> $OUT_FILE
echo -e "\n\n" >> $OUT_FILE


echo -e "bowtie -t -v 3 -a -m 10 --sam -p $NODES $REF_SEQ -1 $READ_1 -2 $READ_2 $OUT_DIR/output_p.sam >> $OUT_FILe\n\n" >> $OUT_FILE
bowtie -t -v 3 -a -m 10 --sam -p $NODES $REF_SEQ -1 $READ_1 -2 $READ_2 $OUT_DIR/output_p.sam >> $OUT_FILE 2>&1

echo -e "\n\nFinished bowtie..." >> $OUT_FILE
date >> $OUT_FILE
echo -e "\n\n" >> $OUT_FILE

CUR_DIR=`pwd`
cd $OUT_DIR
samtools view -h -b -S -o output_p.bam output_p.sam 
samtools sort output_p.bam bowtie-sorted
samtools index bowtie-sorted.bam
rm output_p.sam
rm output_p.bam
cd $CUR_DIR

echo -e "\n\nCompleted BAM file conversion..." >> $OUT_FILE
date >> $OUT_FILE
echo -e "\n\n" >> $OUT_FILE

#########################################################
# Processing of unmapped paired reads:

#bowtie -t -v 3 -a -m 10 --sam -p $NODES $REF_SEQ -1 $READ_1 -2 $READ_2 $OUT_DIR/output_p.sam --un $OUT_DIR/unmapped_p.fq >> $OUT_FILE 2>&1

#echo -e "\n" >> $OUT_FILE
#echo -e "bowtie -t -v 3 -a -m 10 --best --sam -p $NODES $REF_SEQ $OUT_DIR/unmapped_p_1.fq #$OUT_DIR/output_1.sam --un $OUT_DIR/unmapped_1.fq >> $OUT_FILE\n\n" >> $OUT_FILE 2>&1
#bowtie -t -v 3 -a -m 10 --best --sam -p $NODES $REF_SEQ $OUT_DIR/unmapped_p_1.fq $OUT_DIR/output_1.sam --un $OUT_DIR/unmapped_1.fq >> $OUT_FILE 2>&1

#echo -e "\n" >> $OUT_FILE 2>&1
#echo -e "bowtie -t -v 3 -a -m 10 --best --sam -p $NODES $REF_SEQ $OUT_DIR/unmapped_p_2.fq #$OUT_DIR/output_2.sam --un $OUT_DIR/unmapped_2.fq >> $OUT_FILE\n\n" >> $OUT_FILE 2>&1
#bowtie -t -v 3 -a -m 10 --best --sam -p $NODES $REF_SEQ $OUT_DIR/unmapped_p_2.fq $OUT_DIR/output_2.sam --un $OUT_DIR/unmapped_2.fq >> $OUT_FILE 2>&1

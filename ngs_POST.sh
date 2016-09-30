#!/bin/bash

# Copyright (c) 2012-2014, Stephen Fisher and Junhyong Kim, University of
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
# OUTPUT: $SAMPLE/trim/unaligned_1.fq.gz
#
# PAIRED-END READS
# INPUT: $SAMPLE/trim/unaligned_1.fq and $SAMPLE/trim/unaligned_2.fq
# OUTPUT: $SAMPLE/trim/unaligned_1.fq.gz and $SAMPLE/trim/unaligned_2.fq.gz
##########################################################################################

##########################################################################################
# USAGE
##########################################################################################

NGS_USAGE+="Usage: `basename $0` post OPTIONS sampleID    --  clean up RUM and trimmed data\n"

##########################################################################################
# HELP TEXT
##########################################################################################

ngsHelp_POST() {
	echo -e "Usage:\n\t`basename $0` post [-i inputDir] [-g group] sampleID"
	echo -e "Input:\n\tsampleID/INPUTDIR/*.fq"
	echo -e "Output:\n\tsampleID/INPUTDIR/*.fq.gz\n"
	echo -e "Compresses all files that end with 'fq'. For example the unaligned_1.fq file in the trim directory will be compressed with gzip and renamed unaligned_1.fq.gz."
	echo -e "With -g option run chgrp -R [group] on sample."
}

##########################################################################################
# LOCAL VARIABLES WITH DEFAULT VALUES. Using the naming convention to
# make sure these variables don't collide with the other modules.
##########################################################################################

ngsLocal_POST_INP_DIR="trim"
ngsLocal_POST_CHGRP=""

##########################################################################################
# PROCESSING COMMAND LINE ARGUMENTS
# POST args: sampleID
##########################################################################################

ngsArgs_POST() {
	if [ $# -lt 1 ]; then printHelp "POST"; fi

	# getopts doesn't allow for optional arguments so handle them manually
	while true; do
		case $1 in
			-i) ngsLocal_POST_INP_DIR=$2
				shift; shift;
				;;
			-g) ngsLocal_POST_CHGRP=$2
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
# POST command. Compresses trimming data.
##########################################################################################

ngsCmd_POST() {
	prnCmd "# BEGIN: POST PROCESSING"
	
	if [[ $ngsLocal_POST_CHGRP ]]; then
	    prnCmd "chgrp -R $ngsLocal_POST_CHGRP $SAMPLE"
	    chgrp -R $ngsLocal_POST_CHGRP $SAMPLE
	fi

	prnCmd "gzip -f $SAMPLE/$ngsLocal_POST_INP_DIR/*fq"
	if ! $DEBUG; then 
		gzip -f $SAMPLE/$ngsLocal_POST_INP_DIR/*fq
	fi
	
	prnCmd "# FINISHED: POST PROCESSING"
}

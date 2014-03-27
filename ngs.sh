#!/bin/bash -e

# Copyright (c) 2012-2014 Stephen Fisher and Junhyong Kim, University of
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
# Adding new processing step:
# 1. create a shell script containing the following elements ("COMMAND" should be a unique command name")
#    section: SOURCE MODULES HERE
#
# 2. ngsUsage_COMMAND  --  one line command usange and description
#    section: ADD MODULE USAGE HERE
#
# 3. ngsHelp_COMMAND  --  expanded command help
#    section: ADD MODULE HELP HERE
#
# 4. ngsArgs_COMMAND()  --  function for processing command line arguments
#    section: ADD MODULE ARGUMENT FUNCTION HERE
#
# 5. ngsCmd_COMMAND()  --  function that performs command operation
#    section: ADD MODULE COMMAND FUNCTIONS HERE
##########################################################################################

VERSION=2.0.0-alpha

# each module should output a tab-delimited list of file and program
# version information. This file should have two lines, the first line
# being header information and the second line being the versions. The
# prnVersion() command should be used to generate this file. This file
# will live in the respective module subdirectory (ie
# $SAMPLE/MODULE/VERSION)
VERSION_FILE="versions"

DEBUG=false   # disable commands when true, use to see what commands would be run.

# output every line of code to the console when running
#if "$DEBUG"; then set -x

# this is the location of the demultiplexed files, with each sample in a separate subdirectory named with the sample ID.
RAW=raw

# this is the place where analyzed data will be stored. Each sample will be put into a separate subdirectory
ANALYZED=analyzed

# location of genomic databases and library files
REPO_LOCATION=/lab/repo/resources
BOWTIE_REPO=$REPO_LOCATION/bowtie
RUM_REPO=$REPO_LOCATION/rum2
STAR_REPO=$REPO_LOCATION/star
HTSEQ_REPO=$REPO_LOCATION/htseq
SNP_REPO=$REPO_LOCATION/snp

MODULES=( "HELP" "INIT" "FASTQC" "BLAST" "BOWTIE" "TRIM" "STAR" "RUM" "RUMSTATUS" "POST" "BLASTDB" "HTSEQ" "SNP" "SPAdes" "RSYNC" "STATS" "PIPELINE" "VERSION" )

# when modules are loaded, the add their usage to this variable.
NGS_USAGE=""

# default is paired-end.
SE=false  

# make comparisons case insensitive. Need to use [[ and ]] in if
# conditionals, rather than [ and ].
shopt -s nocasematch

# cause an error to happen if trying to use an unset variable. We
# don't use the -u option when launching bash above as that may cause
# bash initialization scripts to error.
set -o nounset

# get OS name. 
OS_VERSION=$(uname)
# If OS isn't "Darwin" (Mac) then assume "Linux" (RedHat /
# Centos). Other OS versions could be added here. 
case ${OS_VERSION} in
	Darwin)
		# Grep on the Mac (and likely BSD) does not have a "-P" option
		# ("perl-regexp"). We use the "-P" option on Linux at various
		# places.
		GREPP="grep"
		;;
	*)
		GREPP="grep -P"
		;;
esac

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# LOAD MODULES
for module in "${MODULES[@]}"; do
	source ngs_${module}.sh
done
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###############################################################################################
# HELPER FUNCTIONS

# print the help text for a module.
printHelp() {
	local NOT_FOUND=true

	for module in "${MODULES[@]}"; do
		if [[ "$module" == "$1" ]]; then
			ngsHelp_$1
			NOT_FOUND=false
			break
		fi
	done

	if $NOT_FOUND; then echo -e $NGS_USAGE; fi

	exit 1
}

# output version information to the module subdirectory. We assume the
# output directory already exists ($SAMPLE/MODULE). We expect three
# arguments: module, header, values.
prnVersion() {
	if [[ $# -ne 3 ]]; then prnError "prnVersion() requires 3 arguments. Only received $#. arguments"; fi

	# we can't rely on COMMAND to know the module calling this
	# function since COMMAND might be pipeline.
	outFile="$SAMPLE/$1/$VERSION_FILE"

	# append pipeline version number
	header="pipeline\t$2"
	values="$VERSION\t$3"

	# we intentionally write over the previous file, if it exists
	echo -e $header > $outFile
	echo -e $values >> $outFile
}

# this will uniformily format the output that is put into the JOURNAL
# file. The following bash command can be used to strip off the time
# stamp and generate a executable bash file:
# cat analysis.log | awk -F\\t '{print \$2}'
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

		echo -ne `date` >> $JOURNAL
		echo -ne "\t" >> $JOURNAL
		echo "##################################################################" >> $JOURNAL
	fi
	echo -ne `date` >> $JOURNAL
	echo -ne "\t" >> $JOURNAL
	echo -ne $1 >> $JOURNAL
	echo >> $JOURNAL
	if [[ $1 == *"# FINISHED"* ]]; then
		echo -ne `date` >> $JOURNAL
		echo -ne "\t" >> $JOURNAL
		echo "##################################################################" >> $JOURNAL

		# insert extra line between sections
		echo >> $JOURNAL

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

# print warning to console
prnWarning() {
	echo -e "\n************************************************" >& 2
	echo -ne "WARNING: " >& 2
	echo -e `date` >& 2
	echo -e $1 >& 2
}

# exit on error
prnError() {
	echo -e "\n************************************************" >& 2
	echo -ne "ERROR: " >& 2
	echo -e `date` >& 2
	echo -e $1 >& 2
	exit 1
}

###############################################################################################
# PROCESS COMMAND ARGUMENTS
 
# if no args then print out usage
if [ $# -lt 1 ]; then
	echo -e $NGS_USAGE
	exit 0
fi

# get command as upper case
COMMAND=$( echo $1 | tr "[a-z]" "[A-Z]" )
shift  # shift removes $1 (ie COMMAND) from the argument list

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# PROCESS MODULE'S ARGUMENT FUNCTION
ngsArgs_${COMMAND} $@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# if we've gotten to this point and $SAMPLE is not set, then something
# went wrong and abort
if [[ -z "$SAMPLE" ]]; then
	echo -e "\n************************************************"
	echo -ne "ERROR: "
	echo -e `date`
	echo "Error processing command arguments."
	echo -e $NGS_USAGE
	exit 1
fi

# remove trailing "/" from $SAMPLE, if present
SAMPLE="${SAMPLE%/}"

# adjust the name of the VERSION file, now that we know the SAMPLE
VERSION_FILE="$SAMPLE.${VERSION_FILE}"

# create output directory
if [[ ! -d $SAMPLE ]]; then
	mkdir $SAMPLE
fi

# create log directory. This needs to happen prior to using the
# prnCmd() function, so we can create the output file ($JOURNAL) that
# is used by prnCmd. This happens even during debugging because this
# directory is where the $JOURNAL file is located by default.
if [[ ! -d $SAMPLE/log ]]; then
	mkdir $SAMPLE/log
fi

# set the location of the JOURNAL file, now that we know the SAMPLE
# and COMMAND. We use a date-time stamp plus the module name. The file
# will be located in the $SAMPLE/log directory. We don't use a ":" in
# the timestamp because Mac file systems don't allow colons in file
# names.
JOURNAL="$SAMPLE/log/$(date +%Y-%m-%d_%H-%M).$COMMAND.log"

# If journal file already exists, then regenerate the filename using
# seconds.
if [ -f $JOURNAL ]; then
	JOURNAL="$SAMPLE/log/$(date +%Y-%m-%d_%H-%M-%S).$COMMAND.log"
fi

# STATS shouldn't write anything to the log file
if [[ "$COMMAND" != "stats" ]]; then
	# log version and run-time information
	if $DEBUG; then prnCmd "# DEBUG MODE"; fi
	_cmd=`basename $0`
	_args=`echo $@`
	prnCmd "# COMMAND: $_cmd $COMMAND $_args"
fi

###############################################################################################
# RUN COMMANDS

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# RUN MODULE'S COMMAND FUNCTION
ngsCmd_${COMMAND}
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


ngs: Next Generation Sequencing Pipeline
========================================

ngs is designed to automate the processing of next generation sequencing data. The master script ngs.sh is used to run the various "modules" that perform individual tasks. Each module should perform a single task. The module name (uppercase) is the command name (lowercase). So for example the 'INIT' module contains the code to perform the 'init' command. 

The TEMPLATE module is an empty, unused module that can be used as an example to build new modules.

Current modules:

		HELP: expanded command help (ex: ngs.sh help blast)
		INIT: prepare read file(s) for processing
		FASTQC: run FastQC
		BLAST: run blast on randomly sampled subset of reads
		BOWTIE: run bowtie on untrimmed reads
		TRIM: trim adapter and poly-A/T contamination
		RUMALIGN: run RUM on trimmed reads
		RUMSTATUS: get status of RUM run
		POST: clean up RUM and trimmed data
		BLASTDB: create blast database from reads
		HTSEQ: run HTSeq on unique mappers from RUM
		RSYNC: copy data to analyzed directory
		STATS: print stats from blast, trimming and RUM
		PIPELINE: run full pipeline
		VERSION: print version information

Modules are designed to be run independently although many modules depend on the output from other modules. For example the TRIM module uses the read files from the 'orig' subdirectory (created by the INIT module) to create the trimmed reads files that it places in the 'trimAT' subdirectory. The RUMALIGN module expects to find reads in the 'trimAT' directory and hence is expected to run after TRIM.

The PIPELINE module is effectively a meta-module and includes the following modules: init, fastqc, blast, trim, rumalign, post, blastdb, htseq, rsync

The expanded help for each module documents the input files, output files, and required programs needed for that module to function. To view the expanded help for a module use the HELP command. For example to get help on the RUMALIGN module: ngs.sh help rumalign


Ancillary
=================

This directory contains files that are not required by any modules provided herein although thay may be helpful to prepare files for use with the pipeline and/or process files after pipeline processing.


Notes
=================

  - Only tested on Linux OS. Will likely work on a Mac. May work on Windows with Cygwin.


To Do
=================

  - GitHub missing trimming scripts. Need to rewrite existing trimming programs and upload to github. Trimming also needs to log version number.
  - Better documentation on how to use the pipeline and how the files relate to one another.
  - Need to untangle POST and perhap create POSTTRIM, POSTBOWTIE and POSTRUM.
  - Deal with rerunning commands. Example if we run INIT twice, then what happens?
  - Allow PIPELINE to pick up where it left off. Example, if INIT, FASTQC, and BLAST directories exist, then should it start with the TRIMMING command?
  - untangle fastqc from trimming. have 2 fastqc commands, one for orig/* and one for trimAT/*.
  - save rum output in rum.$SPECIES instead of rum.trim. This would have implications in rumalign, rumstatus, htseq, stats, and post
  - add flag to STATS that will run stats on all samples and output stats to the specified xls file 
  - need updated trimming scripts. They need to offer version information
  - make journal output optional
  - add argument for setting of journal file
  - make DEBUG an argument
  - Need to version each module.
  - Integrate annotateGeneCnts.py into HTSEQ
  - Update STATS to include percent of unique reads that overlap


No Warranty
=================

Unless otherwise noted in the individual applications, the following disclaimer applies to all applications provided herein.

There is no warranty to the extent permitted by applicable law. Except when otherwise stated in writing the copyright holders and/or other parties provide these applications "as is" without warranty of any kind, either expressed or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. The entire risk as to the quality and performance of these applications, and data is with you. Should these applications or data prove defective, you assume the cost of all necessary servicing, repair or correction.

In no event unless required by applicable law or agreed to in writing will any copyright holder, or any other party who may modify and/or redistribute these applications as permitted above, be liable to you for damages, including any general, special, incidental or consequential damages arising out of the use or inability to use these applications and data (including but not limited to loss of data or data being rendered inaccurate or losses sustained by you or third parties or a failure of these applications and data to operate with any other programs), even if such holder or other party has been advised of the possibility of such damages.

ngs: Next Generation Sequencing Pipeline
========================================

ngs is designed to automate the processing of next generation sequencing data. The master script ngs.sh is used to run the various "modules" that perform individual tasks. Each module should perform a single task. The module name (uppercase) is the command name (lowercase). So for example the 'INIT' module contains the code to perform the 'init' command. 

The TEMPLATE module is an empty, unused module that can be used as an example to build new modules.

Current modules:
		HELP
		INIT
		FASTQC
		BLAST
		BOWTIE
		TRIM
		RUMALIGN
		RUMSTATUS
		POST
		BLASTDB
		HTSEQ
		RSYNC
		STATS
		PIPELINE
		VERSION


* Requires these external programs:
  - RUM
  - samtools
  - blast
  - bowtie
  - fastqc
  - various Linux programs such as grep and rsync
  
* Only tested on Linux OS. Will likely work on a Mac. Will probably not work on Windows.

* ToDo:
  - GitHub missing trimming scripts. Need to rewrite existing trimming programs and upload to github.
  - Better documentation on how to use the pipeline and how the files relate to one another.
  - Need to untangle POST and perhap create POSTTRIM and POSTRUM.


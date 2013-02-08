ngs: Next Generation Sequencing Pipeline
========================================

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
  - Rewrite SOP.sh to be more module, better handling of arguments, more clear argument options, better handling of single versus pair-end reads.

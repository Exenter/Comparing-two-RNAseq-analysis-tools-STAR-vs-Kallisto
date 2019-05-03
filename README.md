# Comparing two RNAseq analysis tools

This repository is the final deliverable for our tutored project as a group of Master bioinformatics students.
Goal was to make and benchmark 2 pipeline using Star or Kallisto for RNAseq analysis. These were made with Snakemake.

# Kallisto pipeline
Pretty straightforward,it runs the 2 main Kallisto commands, the first one indexes the fastas, the second runs the quantification tool.

# Star pipeline
A bit more sofisticated.
Here we have to use 3 different piece of software :
  - Star itself, to index the genome and map the reads on it
  - Samtools, to convert the .sam ouput of Star into a .bam, and then to sort and index this file
  - Htseq-count, which is the quantification tool for the reads we got using Star.

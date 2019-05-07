################################
######## Pipeline STAR #########
################################
# Les fichiers de génomes doivent être dans /genome/
# L'extension doit être .fa.gz
# Les reads doivent être dans /reads/
# Leur extension est .fastq.gz
#
# Pour des données paired-end (2 fichiers : _read1/_read2.fastq):
# input:
#     genome.fa.gz
#     _read1.fastq
#     _read2.fastq
# output:
#     ?
#
# processus :
#     Index du génome de référence
#     Mapping de données (output.sam)
#     sam to bam
#     Rangement des données par coordonnées génomique (samtools sort)
#     Indexer le fichier ordonnée (samtools index)
#     Comptage des reads in feature (htseq_cout)
#
################################
################################
PRE = "HS_chr22"
GENOME = "input/{0}.fa".format(PRE)
FAIDX = "faidx/{0}.fai".format(PRE)
GTF = "input/hs_grch37.gtf".format(PRE)
INDEX = "genome/SAindex".format(PRE)
BAM = "bam/Aligned.out.bam".format(PRE)
BAM_SORTED = "bam/Aligned.out.sorted.bam".format(PRE)
BAM_SORTED_INDEX = "bam/Aligned.out.sorted.bam.bai".format(PRE)

# Fichiers de sortie
rule all:
  input:
    "out/counts_htseq.txt"
   
#expand(out/counts{sample}.txt, sample = qqch (on peut aller chercher dans config))
#Sample = C1, C2 pour la liste des préfixes, à mettre en haut ou dans config
#"Samples" : [ "C1", "C2"]
#expand("{prefix}.count.txt", prefix = config.samples)

# index du génome de référence avec STAR
rule star_index:
  input:
    GENOME = GENOME
  output:
    INDEX
  benchmark:
    repeat("benchmarks/benchmark.star_index.txt",3)
  message:
    "Indexing reference genome ..."
  shell:
      "STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ./genome --genomeFastaFiles ./{input}"
    


 #Mapping
rule map_reads:
  input:
    INDEX = INDEX,
    R1 = "input/reads_1.fastq", #{prefix}_1.fastq => préciser dans config que le séparateur est _ etc
    R2 = "input/reads_2.fastq"
  output:
    "Aligned.out.sam"
  benchmark:
    repeat("benchmarks/benchmark.star_map_reads.txt",3)
  message:
    "Mapping reads on the genome ..."
  shell:
    "STAR --genomeDir ./genome --readFilesIn ./{input.R1} ./{input.R2}"

## Samtools faidx
rule samtools_faidx:
  input:
    GENOME = GENOME
  output:
    FAIDX
  benchmark:
    repeat("benchmarks/benchmark.star_samtools_faidx.txt",3)
  message:
    "Indexing..."
  shell:
    "samtools faidx {input} > {output}"


## Conversion .sam to .bam
rule samTobam:
  input:
    SAM = "Aligned.out.sam",
    FAIDX = FAIDX
  output:
    BAM
  benchmark:
    repeat("benchmarks/benchmark.samtools_samtobam.txt",3)
  message:
    "Converting .sam to .bam ..."
  shell:
    "samtools view -bt {input.FAIDX} -S {input.SAM} > {output}"

# #
## Samtools sort
rule samtools_sort:
  input:
    BAM = BAM
  output:
    BAM_SORTED
  benchmark:
    repeat("benchmarks/benchmark.star_samtools_sort.txt",3)
  message:
    "Sorting bam file ..."
  shell:
    "samtools sort {input} -o  {output[0]}"

# # #
## Samtools Index
rule samtools_index:
  input:
    BAM_SORTED
  output:
    BAM_SORTED_INDEX
  benchmark:
    repeat("benchmarks/benchmark.star_samtools_index.txt",3)
  message:
    " Indexing bam file ..."
  shell:
    "samtools index {input}"
#
## Comptage des reads
rule htseq_count:
  input:
    BAM_SORTED = BAM_SORTED,
    GTF = GTF
  output:
    "out/counts_htseq.txt" #htseq counts et kallisto counts
  benchmark:
    repeat("benchmarks/benchmark.star_htseqcount.txt",3)
  message:
    "Counting reads"
  shell:
    "htseq-count -f bam -s reverse -r pos {input.BAM_SORTED} {input.GTF} > {output}"

#organiser dossiers et wildcards

#utiliser des fichiers .rules et une fichier snakefile qui a juste les deux rules
#fichiers de config: json, dictionnaire python avec clés etc, en haut du fichier
#snakefile il faut mettre config et il va vérifier dans le dossier s'il n'existe pas 
#un fichier config, qui aura les noms de fichier genome etc
#au lieu d'aller chercher le GTF, on ira chercher le config.GTF avec dans le fichier
#config {"GTF" : "hg19.gtf"}
#argument benchmark dans chaque règle en donnant fichier de sortie
#on peut lancer 10 fois chaque règle pour faire du benchmarking, puis faire les stats dessus

#barplot kallisto vs STAR + htseq count => sur le barplot faire temps cumulé sur chaque règle avec ggplot2

#-s reverse pour htseq count selon fichiers (truseq)
#on peut mettre no par défaut dans conf car c'est le plus fréquent

#mettre sur github avec read me pour livrable: ça sera le livrable
#read me détaille chaque script comment ça se lance, dans le fichier config quoi correspond
#à quoi, commenter les fichiers

#trouver option sur kallisto pour floor
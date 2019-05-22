###########################################################
######## Pipeline STAR - STAR QUANTMODE - KALLISTO#########
###########################################################
# Put the Snakefile in a folder along with the config.json, and an /input/ folder
# Reads must be in the /input/ folder
# You must change the read extension in the config file, ex: "READ_EXTENSION" : ".fastq"
# You must change the paths for the DNA, cDNA and the GTF file in the config.json file, ending with a /
# The genomes must be .fa
# If the reads are compressed, change the SUP_STAR_COMPR parameter in the config file to "--readFilesCommand zcat"
# You must add the list of experimental conditions as a list in the config file, ex: "SAMPLE": ["C1reads", "C2reads"]
# You can change the separator between the read name and number, ex: "C1reads_1" "C1reads_2"
#
# For paired-end reads (2 files)
#
# To launch: move to your directory and type snakemake in the terminal, add -j [number of threads] to choose the total number of threads
#
# Process :
#     Index the reference genome with STAR
#     Map the reads (output: sam) with STAR
#     sam to bam
#     Sort by position (samtools sort)
#     Index the sorted file (samtools index)
#     Count the reads by feature (htseq_count) or with STAR quantmode
#
#     Index the cDNA genome with Kallisto
#     Quantification with Kallisto
#
#
################################
################################
configfile: "config.json"

PRE_DNA = config["DNA"]
PRE_CDNA = config["CDNA"]
PRE_GTF = config["GTF"]

DNA_PATH = config["DNA_PATH"]
CDNA_PATH = config["CDNA_PATH"]
GTF_PATH = config["GTF_PATH"]

FASTA = DNA_PATH+"/{0}.fa".format(PRE_DNA)
CDNA = CDNA_PATH+"/{0}.fa".format(PRE_CDNA)
GTF = GTF_PATH+"/{0}.gtf".format(PRE_GTF)

FAIDX = DNA_PATH+"/{0}.fai".format(PRE_DNA)
INDEXS = DNA_PATH+"/chrName.txt"
INDEXK = DNA_PATH+"/{0}.idx".format(PRE_CDNA)

BAM = "bam/Aligned.out.bam"
BAM_SORTED = "bam/Aligned.out.sorted.bam"
BAM_SORTED_INDEX = "bam/Aligned.out.sorted.bam.bai"

SEP = config["SEPARATOR"]
R1 = "1"
R2 = "2"
COMPR = config["SUP_STAR_COMPR"]
READ_EXT = config["READ_EXTENSION"]

THREADS = config["THREADS"]

# Output files
rule all:
	input:
		expand("star/out/{sample}counts_htseq.txt", sample = config["SAMPLE"]),
		expand("kallisto/out/{sample}/abundance_genes.tsv", sample = config["SAMPLE"]),
		expand("star/quantmode/{sample}ReadsPerGene.out.tab", sample = config["SAMPLE"])

#index du génome de référence avec STAR
rule star_index:
	input:
		FASTA = FASTA,
		GTF = GTF
	output:
		INDEXS
	params:
		overhang = config["READ_LENGTH_MINUS_1"],
		GENDIR=DNA_PATH
	threads: THREADS
	benchmark:
		repeat("benchmarks/benchmark.star_index.txt",1)
	message:
		"Indexing reference genome ..."
	shell:
		"STAR --runThreadN {threads} --runMode genomeGenerate --sjdbGTFfile {input.GTF} --sjdbOverhang {params.overhang} --genomeDir {params.GENDIR} --genomeFastaFiles {input.FASTA}"
	
 #Mapping
rule star_map_reads:
	input:
		INDEX = INDEXS,
		R1 = "input/{sample}"+SEP+R1+READ_EXT,
		R2 = "input/{sample}"+SEP+R2+READ_EXT
	output:
		"star/sam/{sample}Aligned.out.sam"
	params:
		PREFIX = "star/sam/{sample}",
		COMPR = COMPR,
		GENDIR=DNA_PATH
	threads: THREADS
	benchmark:
		repeat("benchmarks/benchmark.star_map_reads_{sample}.txt",10)
	message:
		"Mapping reads on the genome ..."
	shell:
		"STAR --genomeDir {params.GENDIR} {params.COMPR} --runThreadN {threads} --outFileNamePrefix {params.PREFIX} --readFilesIn ./{input.R1} ./{input.R2}"

#STAR quantmode
rule star_quantmode:
	input:
		INDEX = INDEXS,
		R1 = "input/{sample}"+SEP+R1+READ_EXT,
		R2 = "input/{sample}"+SEP+R2+READ_EXT
	output:
		"star/quantmode/{sample}ReadsPerGene.out.tab"
	params:
		PREFIX = "star/quantmode/{sample}",
		COMPR = COMPR,
		GTF = GTF,
		GENDIR=DNA_PATH
	threads: THREADS
	benchmark:
		repeat("benchmarks/benchmark.star_quantmode{sample}.txt",10)
	message:
		"Mapping and couting reads on the genome ..."
	shell:
		"STAR --genomeDir {params.GENDIR} {params.COMPR} --sjdbGTFfile {params.GTF} --quantMode GeneCounts --runThreadN {threads} --outFileNamePrefix {params.PREFIX} --readFilesIn ./{input.R1} ./{input.R2}"


## Samtools faidx
rule samtools_faidx:
	input:
		FASTA = FASTA
	output:
		FAIDX
	benchmark:
		repeat("benchmarks/benchmark.star_samtools_faidx.txt",1)
	message:
		"Indexing..."
	shell:
		"samtools faidx {input} > {output}"


## Conversion .sam to .bam
rule samTobam:
	input:
		SAM = "star/sam/{sample}Aligned.out.sam",
		FAIDX = FAIDX
	output:
		"star/bam/{sample}Aligned.out.bam"
	benchmark:
		repeat("benchmarks/benchmark.samtools_samtobam_{sample}.txt",10)
	threads: THREADS
	message:
		"Converting .sam to .bam ..."
	shell:
		"samtools view -@ {threads} -bt {input.FAIDX} -S {input.SAM} > {output}"

##
## Samtools sort
rule samtools_sort:
	input:
		"star/bam/{sample}Aligned.out.bam"
	output:
		"star/bam/{sample}Aligned.out.sorted.bam"
	benchmark:
		repeat("benchmarks/benchmark.star_samtools_sort_{sample}.txt",10)
	threads: THREADS
	message:
		"Sorting bam file ..."
	shell:
		"samtools sort -@ {threads} {input} -o  {output}"

# # #
## Samtools Index
rule samtools_index:
	input:
		"star/bam/{sample}Aligned.out.sorted.bam"
	output:
		"star/bam/{sample}Aligned.out.sorted.bam.bai"
	benchmark:
		repeat("benchmarks/benchmark.star_samtools_index_{sample}.txt",10)
	threads: THREADS
	message:
		" Indexing bam file ..."
	shell:
		"samtools index -@ {threads} {input}"

## Read counting
rule htseq_count:
	input:
		BAM_SORTED = "star/bam/{sample}Aligned.out.sorted.bam",
		GTF = GTF
	output:
		"star/out/{sample}counts_htseq.txt" 
	benchmark:
		repeat("benchmarks/benchmark.star_htseqcount_{sample}.txt",5)
	message:
		"Counting reads"
	shell:
		"htseq-count -f bam -s reverse -r pos {input.BAM_SORTED} {input.GTF} > {output}"

##kallisto

rule kallisto_index:
	input:
		CDNA = CDNA
	output:
		INDEXK
	benchmark:
		repeat("benchmarks/benchmark.kallisto_index.txt",1)
	shell:
		"kallisto index -i {output} {input}"

rule kallisto_quant:
	input:
		R1 = "input/{sample}"+SEP+R1+READ_EXT,
		R2 = "input/{sample}"+SEP+R2+READ_EXT,
		INDEX = INDEXK
	params:
		outdire="kallisto/out/{sample}"
	threads: THREADS
	output:
		"kallisto/out/{sample}/abundance.h5",
		"kallisto/out/{sample}/abundance.tsv",
		"kallisto/out/{sample}/run_info.json"
	benchmark:
		repeat("benchmarks/benchmark.kallisto_quant_{sample}.txt",10)
	shell:
		"kallisto quant "
		"-t {threads} "
		"-i {input.INDEX} "
		"-b 30 "
		"-o {params.outdire} "
		"{input.R1} {input.R2}"
## converting transcript counts to gene counts
rule quant_to_gene:
	input:"kallisto/out/{sample}/abundance.h5"
	output:"kallisto/out/{sample}/abundance_genes.tsv"
	params:
		IDs = config["list_gene_transcripts"],
		col = "Transcript stable ID version"
	benchmark:
		repeat("benchmarks/benchmark.quant_to_gene_{sample}.txt",10)
	script:
		config["TRANS_TO_GENE_CONVERSION"]

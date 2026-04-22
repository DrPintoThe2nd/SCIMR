import os

#mamba env create -f env_scimr.yml
#mamba create -n scimr minimap2 freebayes seqkit seqtk samtools sambamba bamtools bcftools vcftools pixy snakemake mosdepth fastqc multiqc bbmap trim-galore

configfile: "config_scimr_v0.2.json"

R1 = config["R1_suffix"]
R2 = config["R2_suffix"]
males = config["males"]
females = config["females"]
samples = males + females
genome = config["genome"]

threads = config["per_job_threads"]
clump_mem = config["clump_mem"]
window = config["window_size"]
ChrNum = config["ChrNum"]

rule all:
	input:
##setup links rule(s)
		expand("fastqs/{sample}.R1.fq.gz", sample=samples),
		expand("fastqs/{sample}.R2.fq.gz", sample=samples),
##fastqc_analysis rule
		expand("fastqs/{sample}.R1_fastqc.html", sample=samples),
		expand("fastqs/{sample}.R2_fastqc.html", sample=samples),
##trim_galore_pe rule
		expand("trimmed/{sample}.R1_val_1_fastqc.html", sample=samples),
		expand("trimmed/{sample}.R2_val_2_fastqc.html", sample=samples),
##remove_PCR_dups rule
		expand("cleaned/{sample}.R1.fq.gz", sample=samples),
		expand("cleaned/{sample}.R2.fq.gz", sample=samples),
##minimap2_bam rule
		expand("mapped_reads/{sample}.bam", sample=samples),
		expand("mapped_reads/{sample}.bam.bai", sample=samples),
		expand("mapped_reads/{sample}.bam.stat", sample=samples),
##mosdepth rule
		expand("cov_F/F_{female}.regions.bed.gz", female=females),
		expand("cov_F/F_{female}.mosdepth.summary.txt", female=females),
		expand("cov_M/{male}.regions.bed.gz", male=males),
		expand("cov_M/{male}.mosdepth.summary.txt", male=males),
##freebayes_setup rule
		"called_vcf/bams.txt",
##freebayes rule
		"called_vcf/freebayes_calls.vcf.gz",
		"called_vcf/freebayes_calls.vcf.gz.tbi",
##multiqc rule
		"multiqc_report.html",
##vcftools prep rule(s)
		"called_vcf/females.txt",
		"called_vcf/males.txt",
##vcftools Fst rule
		"SCIMR_M_F.windowed.weir.fst",
##plot_Fst rule
        "SCIMR_M_F.windowed.weir.fst.png",
        "SCIMR_M_F.windowed.weir.fst.pdf",

rule link_fastq_R1:
	input:
		"{sample}" + R1
	output:
		"fastqs/{sample}.R1.fq.gz"
	shell:
		"""
		mkdir -p fastqs
		ln -sf $(realpath {input}) {output}
		sleep 1
		"""

rule link_fastq_R2:
	input:
		"{sample}" + R2
	output:
		"fastqs/{sample}.R2.fq.gz",
	shell:
		"""
		mkdir -p fastqs
		ln -sf $(realpath {input}) {output}
		sleep 1
		"""

rule fastqc_analysis:
	input:
		r1 = "fastqs/{sample}.R1.fq.gz",
		r2 = "fastqs/{sample}.R2.fq.gz",
	output:
		out1 = "fastqs/{sample}.R1_fastqc.html",
		out2 = "fastqs/{sample}.R2_fastqc.html"
	threads: threads
	shell:
		"""
		fastqc -t {threads} {input.r1} {input.r2} -o fastqs/;
		"""

rule trim_galore_pe:
	input:
		r1 = "fastqs/{sample}.R1.fq.gz",
		r2 = "fastqs/{sample}.R2.fq.gz",
	output:
		r1_trim = "trimmed/{sample}.R1_val_1.fq.gz",
		r2_trim = "trimmed/{sample}.R2_val_2.fq.gz",
		html1   = "trimmed/{sample}.R1_val_1_fastqc.html",
		html2   = "trimmed/{sample}.R2_val_2_fastqc.html",
	params:
		outdir = "trimmed",
	threads: threads
	shell:
		"""
		mkdir -p {params.outdir}
		trim_galore --paired --cores {threads} --fastqc {input.r1} {input.r2} -o {params.outdir}
		"""

rule remove_PCR_dups:
	input:
		r1 = "trimmed/{sample}.R1_val_1.fq.gz",
		r2 = "trimmed/{sample}.R2_val_2.fq.gz",
	output:
		out1 = "cleaned/{sample}.R1.fq.gz",
		out2 = "cleaned/{sample}.R2.fq.gz",
	params:
		outdir = "cleaned",
		mem = clump_mem,
	shell:
		"""
		mkdir -p {params.outdir}
		clumpify.sh -Xmx{params.mem}g dedupe=t in={input.r1} in2={input.r2} out={output.out1} out2={output.out2} tmpdir={params.outdir} -da showspeed=t overwrite=t 2>/dev/null
		"""

rule minimap2_bam:
	input:
		r1 = "cleaned/{sample}.R1.fq.gz",
		r2 = "cleaned/{sample}.R2.fq.gz",
	output:
		bam  = "mapped_reads/{sample}.bam",
		bai  = "mapped_reads/{sample}.bam.bai",
		stat = "mapped_reads/{sample}.bam.stat",
	params:
		outdir = "mapped_reads",
		header = str("-R '@RG\\tID:{sample}\\tSM:{sample}'"),
		genome = genome,
	threads: threads
	shell:
		"""
		mkdir -p {params.outdir}
		minimap2 -ax sr -t{threads} {params.header} {params.genome} {input.r1} {input.r2} | samtools fixmate -@ {threads} - - | samtools sort -@ {threads} -O bam - -o {output.bam};
		sambamba index -t{threads} {output.bam};
		sambamba flagstat -t{threads} {output.bam} > {output.stat};
		"""

rule calc_cov_F:
	input:
		bam  = "mapped_reads/{female}.bam",
	output:
		cov = "cov_F/F_{female}.regions.bed.gz",
		sum = "cov_F/F_{female}.mosdepth.summary.txt",
	params:
		prefix = lambda wc: "cov_F/F_" + wc.female,
		window = window,
	threads: threads
	shell:
		"""
		mkdir -p cov_F
		MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{threads} {params.prefix} {input.bam}
		"""

rule calc_cov_M:
	input:
		bam  = "mapped_reads/{male}.bam",
	output:
		cov = "cov_M/{male}.regions.bed.gz",
		sum = "cov_M/{male}.mosdepth.summary.txt",
	params:
		prefix = lambda wc: "cov_M/" + wc.male,
		window = window,
	threads: threads
	shell:
		"""
		mkdir -p cov_M
		MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{threads} {params.prefix} {input.bam}
		"""

rule freebayes_setup:
	input:
		bams = expand("mapped_reads/{sample}.bam", sample=samples),
	output:
		"called_vcf/bams.txt",
	params:
		indir = "mapped_reads",
		outdir = "called_vcf",
	shell:
		"""
		mkdir -p {params.outdir};
		ls {params.indir}/*.bam > {output};
		"""

rule freebayes:
	input:
		"called_vcf/bams.txt",
	output:
		vcf = "called_vcf/freebayes_calls.vcf.gz",
		tbi = "called_vcf/freebayes_calls.vcf.gz.tbi",
	params:
		genome = genome,
	shell:
		"""
#		freebayes -f {params.genome} -g 100 --report-monomorphic --gvcf -L {input} | bgzip -c > {output.vcf};
		freebayes -f {params.genome} -g 100 -L {input} | bgzip -c > {output.vcf};
		tabix {output.vcf};
		"""

##runs multiqc for all types of reads used and generated above
rule multiqc_all:
	input:
		expand("trimmed/{sample}.R1_val_1_fastqc.html", sample=samples),
	output:
		"multiqc_report.html",
	shell:
		"""
		multiqc .
		"""

rule Fst_prep:
	input:
		checkpoint = "called_vcf/bams.txt",
	output:
		females = "called_vcf/females.txt",
		males = "called_vcf/males.txt",
	params:
		config_F = females,
		config_M = males,
	shell:
		"""
		echo {params.config_F} | sed 's/[[:space:]]/\\n/g' > {output.females};
		echo {params.config_M} | sed 's/[[:space:]]/\\n/g' > {output.males};
		"""

rule calc_Fst:
	input:
		vcf = "called_vcf/freebayes_calls.vcf.gz",
		tbi = "called_vcf/freebayes_calls.vcf.gz.tbi",
		females = "called_vcf/females.txt",
		males = "called_vcf/males.txt"
	output:
		"SCIMR_M_F.windowed.weir.fst",
	params:
		window = window,
	shell:
		"""
		vcftools --gzvcf {input.vcf} --fst-window-size {params.window} --fst-window-step {params.window} --weir-fst-pop {input.males} --weir-fst-pop {input.females}  --stdout > {output}
		"""

rule plot_Fst:
	input:
		"SCIMR_M_F.windowed.weir.fst",
	output:
		png = "SCIMR_M_F.windowed.weir.fst.png",
		pdf = "SCIMR_M_F.windowed.weir.fst.pdf",
	params:
		chrnum = ChrNum,
	shell:
		"""
        Rscript manhattan_Fst.R {input} {input} {params.chrnum}
		"""

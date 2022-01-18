rule sniffles_call:
	input:
		bam="alignments/{sample}.{aligner}.srt.bam",
		bai="alignments/{sample}.{aligner}.srt.bam.bai"
	output:
		"results/{sample}.{aligner}.sniffles.vcf"
	threads: config["threads"]
	log:
		"logs/{sample}.{aligner}.sniffles.log"
	conda:
		"../envs/sniffles.yaml"
	shell:
		"sniffles -m {input.bam} -v {output} --genotype -s 2 -l 50 -t {threads} 2>{log}"

rule cuteSV_call:
	input:
		ref=config["genome"],
		bam="alignments/{sample}.{aligner}.srt.bam",
		bai="alignments/{sample}.{aligner}.srt.bam.bai"
	output:
		"results/{sample}.{aligner}.cutesv.vcf"
	threads: config["threads"]
	log:
		"logs/{sample}.{aligner}.cutesv.log"
	params:
		name="{sample}",
		wd="results/{sample}_{aligner}_cutesv"
	conda:
		"../envs/cutesv.yaml"
	shell:
		"mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} -S {params.name} {input.bam} {input.ref} {output} {params.wd} 2>{log}"


rule clean_sniffles:
	input:
		exclude=config["exclude"],
		vcf="results/{sample}.{aligner}.sniffles.vcf"
	output:
		"results/{sample}.{aligner}.sniffles.srt.pass.vcf"
	threads: config["threads"]
	log:
		"logs/{sample}.{aligner}.sniffles_srt_log"
	conda:
		"../envs/sniffles.yaml"
	shell:
		"""
		cat <(cat {input.vcf} | grep "^#") <(cat {input.vcf} | grep -vE "^#" | grep -v "0/0:" | sort -k1,1 -k2,2g) |bcftools view -f PASS |grep -w -v -f {input.exclude} > {output}
        """    

rule clean_cutesv:
	input:
		exclude=config["exclude"],
		vcf="results/{sample}.{aligner}.cutesv.vcf"
	output:
		"results/{sample}.{aligner}.cutesv.srt.pass.vcf"
	threads: config["threads"]
	log:
		"logs/{sample}.{aligner}.cutesv_srt_log"
	conda:
		"../envs/cutesv.yaml"
	shell:
		"""
		bcftools view -f PASS {input.vcf} | grep -w -v -f {input.exclude} |bcftools sort -o {output} - 
        """

rule sv_merge_hc:
	input:
		expand(f"results/{{sample}}.{{aligner}}.{{tool}}.srt.pass.vcf", sample=config["samples"].values(),aligner=config["aligners"].values(), tool=config["tools"].values())
	output:
		expand(f"results/somatic.{{aligner}}.vcf", aligner=config["aligners"].values()),
		expand(f"results/{{sample}}.{{aligner}}.hc.vcf", sample=config["samples"].values(), aligner=config["aligners"].values())
	threads: config["threads"]
	log:
		"logs/log.hc.log"
	conda:
		"../envs/survivor.yaml"
	params:
		sample=' '.join(list(config['samples'].values())),
		algn=' '.join(list(config['aligners'].values())),
		germ="results/germ.txt",
		outdir="results"
	shell:
		"""
		for a in {params.algn}; do for s in {params.sample}; do ls {params.outdir}/${{s}}.${{a}}*.srt.pass.vcf > {params.outdir}/${{s}}.${{a}}.list.txt && SURVIVOR merge {params.outdir}/${{s}}.${{a}}.list.txt 1000 2 1 1 0 50 {params.outdir}/${{s}}.${{a}}.tmp.vcf && rm {params.outdir}/${{s}}.${{a}}.list.txt && bcftools view -i 'GT[*]!="./."' -o {params.outdir}/${{s}}.${{a}}.hc.vcf {params.outdir}/${{s}}.${{a}}.tmp.vcf && rm {params.outdir}/${{s}}.${{a}}.tmp.vcf; done && ls {params.outdir}/*.${{a}}.hc.vcf > {params.outdir}/${{a}}.total.list.txt && SURVIVOR merge {params.outdir}/${{a}}.total.list.txt 1000 1 1 1 0 50 {params.outdir}/total.${{a}}.vcf && rm {params.outdir}/${{a}}.total.list.txt && bcftools query -l {params.outdir}/total.${{a}}.vcf |grep G > {params.germ} && bcftools view -e 'GT[@{params.germ}]!="./."' -o {params.outdir}/somatic.${{a}}.vcf {params.outdir}/total.${{a}}.vcf; done 
		"""

rule make_upset:
	input:
		expand(f"results/somatic.{{aligner}}.vcf", aligner=config["aligners"].values())
	output:
		"results/upset.tsv"
	threads: config["threads"]
	log: "logs/upset.log"
	conda:
		"../envs/survivor.yaml"
	params:
		tmptxt="results/tmp.txt",
		mergedvcf="results/merged.vcf"
	shell:
		"""
		bcftools view -h {input} | tail -1 | cut -f 1-3,10- | sed "s/POS/SVTYPE/g" | sed "s/ID/SVLEN/g" > {output}  && \
		bcftools query -f '%CHROM\t%SVTYPE\t%SVLEN[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {input} | awk '{{FS=OFS="\t"}}{{for (i=4; i<=NF; i++) if ($i=="./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $i="0"; else $i="1"}}1' >> {output} 2>{log}
		"""

rule plot_upset:
	input:
		"results/upset.tsv"
	output:
		"results/upset.pdf"
	threads: config["threads"]
	log: "logs/upset_plot.log"
	conda:
		"../envs/plot.yaml"
	params:
		script="workflow/scripts/upset.R"
	shell:
		"Rscript {params.script} {input} {output}"


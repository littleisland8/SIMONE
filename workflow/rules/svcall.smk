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

rule sv_3X:
	input:
		expand(f"results/{{sample}}.{{aligner}}.hc.vcf", sample=config["samples"].values(), aligner=config["aligners"].values())
	output:
		expand(f"results/3X/somatic.{{aligner}}.3x.vcf", aligner=config["aligners"].values())
	threads: config["threads"]
	log:
		"logs/log.hc.3x.log"
	conda:
		"../envs/survivor.yaml"
	params:
		sample=' '.join(list(config['samples'].values())),
		algn=' '.join(list(config['aligners'].values())),
		germ="results/germ.txt",
		hcdir="results",
		outdir="results/3X"
	shell:
		"""
		for a in {params.algn}; do for s in {params.sample}; do bcftools view -i "FORMAT/DR[0:1] >= 3 && FORMAT/DR[1:1] >= 3" -o {params.outdir}/${{s}}.${{a}}.hc.3x.vcf {params.hcdir}/${{s}}.${{a}}.hc.vcf; done && ls {params.outdir}/*.${{a}}.hc.3x.vcf > {params.outdir}/${{a}}.total.list.txt && SURVIVOR merge {params.outdir}/${{a}}.total.list.txt 1000 1 1 1 0 50 {params.outdir}/total.${{a}}.3x.vcf && rm {params.outdir}/${{a}}.total.list.txt && bcftools view -e 'GT[@{params.germ}]!="./."' -o {params.outdir}/somatic.${{a}}.3x.vcf {params.outdir}/total.${{a}}.3x.vcf; done 
		"""

rule sv_5X:
	input:
		expand(f"results/{{sample}}.{{aligner}}.hc.vcf", sample=config["samples"].values(), aligner=config["aligners"].values())
	output:
		expand(f"results/5X/somatic.{{aligner}}.5x.vcf", aligner=config["aligners"].values())
	threads: config["threads"]
	log:
		"logs/log.hc.5x.log"
	conda:
		"../envs/survivor.yaml"
	params:
		sample=' '.join(list(config['samples'].values())),
		algn=' '.join(list(config['aligners'].values())),
		germ="results/germ.txt",
		hcdir="results",
		outdir="results/5X"
	shell:
		"""
		for a in {params.algn}; do for s in {params.sample}; do bcftools view -i "FORMAT/DR[0:1] >= 5 && FORMAT/DR[1:1] >= 5" -o {params.outdir}/${{s}}.${{a}}.hc.5x.vcf {params.hcdir}/${{s}}.${{a}}.hc.vcf; done && ls {params.outdir}/*.${{a}}.hc.5x.vcf > {params.outdir}/${{a}}.total.list.txt && SURVIVOR merge {params.outdir}/${{a}}.total.list.txt 1000 1 1 1 0 50 {params.outdir}/total.${{a}}.5x.vcf && rm {params.outdir}/${{a}}.total.list.txt && bcftools view -e 'GT[@{params.germ}]!="./."' -o {params.outdir}/somatic.${{a}}.5x.vcf {params.outdir}/total.${{a}}.5x.vcf; done 
		"""

rule sv_10X:
	input:
		expand(f"results/{{sample}}.{{aligner}}.hc.vcf", sample=config["samples"].values(), aligner=config["aligners"].values())
	output:
		expand(f"results/10X/somatic.{{aligner}}.10x.vcf", aligner=config["aligners"].values())
	threads: config["threads"]
	log:
		"logs/log.hc.10x.log"
	conda:
		"../envs/survivor.yaml"
	params:
		sample=' '.join(list(config['samples'].values())),
		algn=' '.join(list(config['aligners'].values())),
		germ="results/germ.txt",
		hcdir="results",
		outdir="results/10X"
	shell:
		"""
		for a in {params.algn}; do for s in {params.sample}; do bcftools view -i "FORMAT/DR[0:1] >= 10 && FORMAT/DR[1:1] >= 10" -o {params.outdir}/${{s}}.${{a}}.hc.10x.vcf {params.hcdir}/${{s}}.${{a}}.hc.vcf; done && ls {params.outdir}/*.${{a}}.hc.10x.vcf > {params.outdir}/${{a}}.total.list.txt && SURVIVOR merge {params.outdir}/${{a}}.total.list.txt 1000 1 1 1 0 50 {params.outdir}/total.${{a}}.10x.vcf && rm {params.outdir}/${{a}}.total.list.txt && bcftools view -e 'GT[@{params.germ}]!="./."' -o {params.outdir}/somatic.${{a}}.10x.vcf {params.outdir}/total.${{a}}.10x.vcf; done 
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


rule make_upset_3X:
	input:
		expand(f"results/3X/somatic.{{aligner}}.3x.vcf", aligner=config["aligners"].values())
	output:
		"results/3X/upset.tsv"
	threads: config["threads"]
	log: "logs/upset.3x.log"
	conda:
		"../envs/survivor.yaml"
	params:
		tmptxt="results/3X/tmp.txt",
		mergedvcf="results/3X/merged.vcf"
	shell:
		"""
		bcftools view -h {input} | tail -1 | cut -f 1-3,10- | sed "s/POS/SVTYPE/g" | sed "s/ID/SVLEN/g" > {output}  && \
		bcftools query -f '%CHROM\t%SVTYPE\t%SVLEN[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {input} | awk '{{FS=OFS="\t"}}{{for (i=4; i<=NF; i++) if ($i=="./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $i="0"; else $i="1"}}1' >> {output} 2>{log}
		"""

rule plot_upset_3X:
	input:
		"results/3X/upset.tsv"
	output:
		"results/3X/upset.pdf"
	threads: config["threads"]
	log: "logs/upset_plot.3x.log"
	conda:
		"../envs/plot.yaml"
	params:
		script="workflow/scripts/upset.R"
	shell:
		"Rscript {params.script} {input} {output}"

rule make_upset_5X:
	input:
		expand(f"results/5X/somatic.{{aligner}}.5x.vcf", aligner=config["aligners"].values())
	output:
		"results/5X/upset.tsv"
	threads: config["threads"]
	log: "logs/upset.5x.log"
	conda:
		"../envs/survivor.yaml"
	params:
		tmptxt="results/5X/tmp.txt",
		mergedvcf="results/5X/merged.vcf"
	shell:
		"""
		bcftools view -h {input} | tail -1 | cut -f 1-3,10- | sed "s/POS/SVTYPE/g" | sed "s/ID/SVLEN/g" > {output}  && \
		bcftools query -f '%CHROM\t%SVTYPE\t%SVLEN[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {input} | awk '{{FS=OFS="\t"}}{{for (i=4; i<=NF; i++) if ($i=="./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $i="0"; else $i="1"}}1' >> {output} 2>{log}
		"""

rule plot_upset_5X:
	input:
		"results/5X/upset.tsv"
	output:
		"results/5X/upset.pdf"
	threads: config["threads"]
	log: "logs/upset_plot.5x.log"
	conda:
		"../envs/plot.yaml"
	params:
		script="workflow/scripts/upset.R"
	shell:
		"Rscript {params.script} {input} {output}"

rule make_upset_10X:
	input:
		expand(f"results/10X/somatic.{{aligner}}.10x.vcf", aligner=config["aligners"].values())
	output:
		"results/10X/upset.tsv"
	threads: config["threads"]
	log: "logs/upset.10x.log"
	conda:
		"../envs/survivor.yaml"
	params:
		tmptxt="results/10X/tmp.txt",
		mergedvcf="results/10X/merged.vcf"
	shell:
		"""
		bcftools view -h {input} | tail -1 | cut -f 1-3,10- | sed "s/POS/SVTYPE/g" | sed "s/ID/SVLEN/g" > {output}  && \
		bcftools query -f '%CHROM\t%SVTYPE\t%SVLEN[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {input} | awk '{{FS=OFS="\t"}}{{for (i=4; i<=NF; i++) if ($i=="./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $i="0"; else $i="1"}}1' >> {output} 2>{log}
		"""

rule plot_upset_10X:
	input:
		"results/10X/upset.tsv"
	output:
		"results/10X/upset.pdf"
	threads: config["threads"]
	log: "logs/upset_plot.10x.log"
	conda:
		"../envs/plot.yaml"
	params:
		script="workflow/scripts/upset.R"
	shell:
		"Rscript {params.script} {input} {output}"

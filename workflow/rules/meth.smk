rule index_fast5:
	input:
		"data/{sample}.pass.fastq.gz"
	output:
		"data/{sample}.pass.fastq.gz.index.readdb"
	threads: 1
	conda: 
		"../envs/nanopolish.yaml"
	log: "logs/{sample}_nanopolish_index.log"
	params:
		name="{sample}",
		f5path=config["fast5_data"]
	shell:
		"""
		var=$(find {params.f5path} -name {params.name} -type d -exec readlink -f {{}} +;) && \
		nanopolish index -d ${{var}} {input} 2>{log}
		"""
rule call_meth:
	input:
		fq="data/{sample}.pass.fastq.gz",
		db="data/{sample}.pass.fastq.gz.index.readdb",
		ref=config["genome"],
		bam="alignments/{sample}.minimap2.srt.bam",
		bai="alignments/{sample}.minimap2.srt.bam.bai"
	output:
		"results/{sample}.methcalls.tsv"
	log:
		"logs/{sample}_nanopolish_methylation.log"
	threads: 10
	conda: 
		"../envs/nanopolish.yaml"
	shell:
		"nanopolish call-methylation -t {threads} -r {input.fq} -b {input.bam} -g {input.ref} > {output} 2>{log}"

rule nanopolish_freq:
        input:
                "results/{sample}.methcalls.tsv"
        output:
                "results/{sample}.methfreq.tsv"
        threads: config["threads"]
        conda:
                "../envs/nanopolish.yaml"
        log:
                "logs/{sample}.nanopolish_freq.log"
        params:
                script="workflow/scripts/calculate_methylation_frequency.py"
        shell:
                "python {params.script} {input} > {output} 2>{log}"

rule make_window:
	input:
		"data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
	output:
		"results/intervals.bed"
	threads:  config["threads"]
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/make_interval.log"
	params:
		bin=config["bin"]
	shell:
		"""
		bedtools makewindows -b <(cut -f1,2 {input} | head -24 |awk 'OFS=FS="\t"''{{print $1, "1", $2}}') -w {params.bin} | sortBed > {output} 2>{log}
		"""	

rule parse_meth_freq:
	input:
		intervals="results/intervals.bed",
		methfreq="results/{sample}.methfreq.tsv"
	output:
		"results/{sample}.methfreq.parsed.bed"
	threads: config["threads"]
	conda:
		"../envs/bedtools.yaml"
	log:
		"logs/{sample}.nanopolish_freq.parse.log"
	params:
		name="{sample}"
	shell:
		"""
		tail -n+2 {input.methfreq} | sortBed | awk 'FS=OFS="\t"''{{print $1,$2,$3,$7}}' > results/{params.name}.test.bed && bgzip results/{params.name}.test.bed && tabix results/{params.name}.test.bed.gz && while IFS=$'\t' read -r chrom start end; do tabix results/{params.name}.test.bed.gz $chrom":"$start"-"$end > results/{params.name}.overlap.bed && if [ -s results/{params.name}.overlap.bed ]; then mean=$(cut -f 4 results/{params.name}.overlap.bed | awk '{{ sum += $1; n++ }} END {{ if (n > 0) print sum / n; }}') && echo -e $chrom"\t"$start"\t"$end"\t"$mean"\t"{params.name} >> {output};fi;done<{input.intervals} && rm results/{params.name}.overlap.bed && rm results/{params.name}.test.bed.gz && rm results/{params.name}.test.bed.gz.tbi 2>{log}
		"""

rule pycometh_CpG:
	input:
		meth="results/{sample}.methcalls.tsv",
		ref=config["genome"]
	output:
		bed9="results/{sample}.CpG.bed9",
		tsv="results/{sample}.CpG.tsv"
	threads: config["threads"]
	conda:
		"../envs/pycometh.yaml"
	log: "logs/{sample}_pycoMeth.cpg.log"
	params:
		name="{sample}"
	shell:
		"pycoMeth CpG_Aggregate -i {input.meth} -f {input.ref} -b {output.bed9} -t {output.tsv} -d 5 -s {params.name} -v -p 2>{log}"

rule pycometh_Interval:
	input:
		CpG="results/{sample}.CpG.tsv",
		ref=config["genome"]
	output:
		bed9="results/{sample}.interval.bed9",
		tsv="results/{sample}.interval.tsv"
	threads: config["threads"]
	conda:
		"../envs/pycometh.yaml"
	log: "logs/{sample}_pycoMeth.interval.log"
	params:
		name="{sample}"
	shell:
		"pycoMeth Interval_Aggregate -i {input.CpG} -f {input.ref} -b {output.bed9} -t {output.tsv} -s {params.name} -v -p 2>{log}"

rule pycometh_Meth_Comp:
	input:
		expand(f"results/{{sample}}.interval.tsv", sample=config["samples"].values())
	output:
		directory("results/pycoMeth")
	threads: config["threads"]
	conda:
		"../envs/pycometh.yaml"
	log:"logs/pycoMethcomp.log"
	params:
		ref=config["genome"],
		include=config["chrs"],
		gff3=config["gff3"],
		outdir="results/pycoMeth",
		coups=' '.join(list(config["merged"].values()))
	shell:
		"""
		mkdir -p {params.outdir} && for c in {params.coups}; do g=$(echo ${{c}} |cut -d "_" -f 2) && s=$(echo ${{c}} |cut -d "_" -f 1) && pycoMeth Meth_Comp -s ${{g}} ${{s}} -i results/${{g}}.interval.tsv results/${{s}}.interval.tsv -f {params.ref} -b results/${{g}}_${{s}}.methcomp.bed9 -t results/${{g}}_${{s}}.methcomp.tmp.tsv && head -1 results/${{g}}_${{s}}.methcomp.tmp.tsv > results/header.txt && grep -w -f {params.include} results/${{g}}_${{s}}.methcomp.tmp.tsv > results/${{g}}_${{s}}.tmp.tsv && cat results/header.txt results/${{g}}_${{s}}.tmp.tsv > results/${{g}}_${{s}}.methcomp.tsv && rm results/${{g}}_${{s}}.methcomp.tmp.tsv results/${{g}}_${{s}}.tmp.tsv results/header.txt && pycoMeth Comp_Report -i results/${{g}}_${{s}}.methcomp.tsv -f {params.ref} -g {params.gff3} -o {params.outdir}/${{g}}_${{s}} 2>>{log}; done 
		"""

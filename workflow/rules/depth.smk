rule mosdepth_depth:
	input:
		bam="alignments/{sample}.{aligner}.srt.bam",
		bai="alignments/{sample}.{aligner}.srt.bam.bai"
	output:
		"alignments/{sample}.{aligner}.mosdepth.global.dist.txt"
	log:
		"logs/{sample}.{aligner}_mosdepth.log"
	conda: "../envs/mosdepth.yaml"
	params:
		prefix="alignments/{sample}.{aligner}"
	threads: config["threads"]
	shell:
		"mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"


rule depth_plot:
	input:
		directory(f"alignments")
	output:
		"alignments/cov_plot.pdf"
	log:
		"logs/depth_plot.log"
	conda:
		"../envs/plot.yaml"
	threads: 1
	params:
		script="workflow/scripts/plotcov.R"
	shell:
		"Rscript {params.script} {input} {output}"

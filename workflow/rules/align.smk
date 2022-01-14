rule minimap2_align:
	input:
		config['genome'],
		"data/{sample}.pass.fastq.gz"
	output:
		"alignments/{sample}.minimap2.srt.bam"
	threads: config["threads"]
	log:
		"logs/{sample}.minimap2_alignment.log"
	conda: "../envs/minimap2.yaml"
	params: 
		fastq="data/{sample}.pass.fastq"
	shell:
		"minimap2 --MD -ax map-ont -t {threads} {input} | samtools sort -@ {threads} -o {output} - 2>{log}"

rule index_bam:
	input:
		"alignments/{sample}.{aligner}.srt.bam"
	output:
		"alignments/{sample}.{aligner}.srt.bam.bai"
	threads:config["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools index -@ {threads} {input}"


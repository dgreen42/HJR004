process call_methylation {
	label "HJR004"
	cpus 4
	memory "16GB"
	publishDir "methylation", mode: "copy"
	input:
		val hold
	output: 
		path "somepath"
	script:
	"""
		nanopolish index -d ${fastq_dir} meth_index.fastq
		echo "hi"
	"""
}

workflow methylation_analysis {
	take:
		fastqs
		bams
		ref_genome
	main:
		call_methylation(fastqs, bams)
}

import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

params.fastq = null 
params.ref_annotation = null
params.ref_genome = null
params.sample_sheet = null
params.ndr = 3

process check_tools {
	label "HJR004"
	cpus 1
	output:
		path "tool_versions.txt"
	script:
	"""
	minimap2 --version | sed 's/^/minimap2,/' >> tool_versions.txt
	samtools --version | head -n1 >> tool_versions.txt
	R --version | head -n1 >> tool_versions.txt
	nanopolish --version | head -n1 >> tool_versions.txt
	"""
}

process index_ref {
	label "HJR004"
	cpus 4
	memory "16 GB"
	input:
		path genome
	output:
		path "index.mmi", emit: index
	script:
	"""
		minimap2 -d "index.mmi" ${genome}
	"""
}

process align_and_map {
	label "HJR004"
	cpus 4
	memory "16 GB"
	publishDir "aln", mode: "copy"
	input:
		each fastq
		val idx
		val genome
	output:
		path "${fastq.baseName}_filtered.bam", emit: bampath
	script:
	"""
		minimap2 -ax splice -uf -k14 ${idx} ${fastq} > ${fastq.baseName}_aln.sam
		samtools view -bT ${genome} -o ${fastq.baseName}_out.bam ${fastq.baseName}_aln.sam
		samtools view -q 10 -m 20 -F 2304 -o ${fastq.baseName}_filtered.bam ${fastq.baseName}_out.bam
	"""
}

process count_transcripts {
	label "HJR004"
	cpus 4
	memory "16 GB"
	publishDir "cts", mode: "copy"
	input:
		val genome
		val anno
		val ndr
		path bams
	output:
		path "bambu_out/HJR004_counts_gene.txt"
		path "bambu_out/HJR004_counts_transcript.txt", emit: transcripts
		path "bambu_out/HJR004_CPM_transcript.txt"
		path "bambu_out/HJR004_extended_annotations.gtf"
		path "bambu_out/HJR004_fullLengthCounts_transcript.txt"
		path "bambu_out/HJR004_uniqueCounts_transcript.txt"
	script:
	"""
		bambu_cts.R ${genome} ${anno} ${ndr} ${bams}
	"""
}

process stage_wise_analysis {
	label "HJR004"
	cpus 4
	memory "16GB"
	publishDir "analysis", mode: "copy"
	input:
		val counts
		val anno
		val sample_sheet
	output:
		path "analysis/de_coefficients.csv"
		path "analysis/de_results.csv"
		path "analysis/de_summary.csv"
		path "analysis/dex_adjusted_pval.csv"
		path "analysis/dex_altsplice.csv"
		path "analysis/dex_isoform_proportions.csv"
	script:
	"""
		diff_splice_stageR.R ${counts} ${anno} ${sample_sheet}
	"""
}

process index_methylation {
	label "HJR004"
	cpus 4
	memory "16GB"
	input:
		val fastq_dir
	output:
		path "meth_index.fastq.index"
		path "meth_index.fastq.index.fai"
		path "meth_index.fastq.fastq.index.gzi"
		path "meth_index.fastq.index.readdb"
	script:
	"""
		nanopolish index -d ${fastq_dir} meth_index.fastq
	"""
}

process methylation_analysis {
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
		echo "hi"
	"""
}

workflow {

	error = null

	if (params.ref_annotation) {
		ref_anno = file(params.ref_annotation, type: "file")
		if (!ref_anno.exists()) {
			error = "reference annotation does not exist"
		}
	} else {
		error = "no reference annotation provided"
	}

	if (params.ref_genome) {
		ref_genome = file(params.ref_genome, type: "file")
		if (!ref_genome.exists()) {
			error = "reference genome does not exist"
		}
	} else {
		error = "no reference genome provided"
	}

	if (params.sample_sheet) {
		sample_sheet = file(params.sample_sheet, type: "file")
		if (!sample_sheet.exists()) {
			error = "sample sheet does not exist"
		}
	} else {
		error = "no sample sheet provided"
	}

	if (error) {
		throw new Exception(error)
	}

	check_tools()
	index_ref(ref_genome)
	fq_ch = channel.fromPath(params.fastq + "*.fastq.gz")
	align_and_map(fq_ch, index_ref.out.index, ref_genome)
	bams = align_and_map.out.bampath.collect(x -> file(x, type: "file"))
	count_transcripts(ref_genome, ref_anno, params.ndr, bams)
	stage_wise_analysis(count_transcripts.out.transcripts, ref_anno, sample_sheet)
}

workflow.onComplete {
	println("All done!")
}
workflow.onError {
	println("Well shit..." + error)
}


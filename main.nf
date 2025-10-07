import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

include {methylation_analysis} from "./workflows/methylation_analysis.nf"

params.fastq = null 
params.bams = null
params.ref_annotation = null
params.ref_genome = null
params.sample_sheet = null
params.acronym_list = null
params.ndr = 3
params.m6A_analysis = false

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
		path "${fastq.baseName}_filtered_sorted.bam", emit: bampath
	script:
	"""
		minimap2 -ax splice -uf -k14 ${idx} ${fastq} > ${fastq.baseName}_aln.sam
		samtools view -q 10 -m 20 -F 2304 -bT ${genome} ${fastq.baseName}_aln.sam | samtools sort -o ${fastq.baseName}_filtered_sorted.bam
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
		path "bambu_out/HJR004_extended_annotations.gtf", emit: ext_anno
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
		val acronym_list
	output:
		path "de_coefficients.csv"
		path "de_results.csv"
		path "de_summary.csv"
		path "de_adjusted_results.csv"
		path "de_adjusted_pvalues.csv"
		path "taxa_to_gene_distribution.csv"
		path "dex_adjusted_pval.csv"
		path "dex_altsplice.csv"
		path "dex_isoform_proportions.csv"
		path "dex_isoform_proportions_nod.csv"
		path "dex_isoform_proportions_irt.csv"
		path "dex_isoform_proportions_mrt.csv"
	script:
	"""
		diff_splice_stageR.R ${counts} ${anno} ${sample_sheet} ${acronym_list}
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

	if (params.acronym_list) {
		acronym_list = file(params.acronym_list, type: "file")
		if (!acronym_list.exists()) {
			error = "acronym list does not exist"
		}
	} else {
		error = "no acronym_list provided"
	}

	if (error) {
		throw new Exception(error)
	}

	check_tools()

	if (!params.bams) {
		index_ref(ref_genome)
		fq_ch = channel.fromPath(params.fastq + "*.fastq.gz")
		align_and_map(fq_ch, index_ref.out.index, ref_genome)
		bams = align_and_map.out.bampath.collect(x -> file(x, type: "file"))
	} else {
		bams = file(params.bams + "*.bam")
	}
	
	if (params.m6A_analysis) {
		me = methylation_analysis(fastq, bams, ref_genome)
	}

	count_transcripts(ref_genome, ref_anno, params.ndr, bams)
	stage_wise_analysis(count_transcripts.out.transcripts, count_transcripts.out.ext_anno, sample_sheet, acronym_list)
}

workflow.onComplete {
	println("All done!")
}
workflow.onError {
	println("Well shit... something went wrong")
}


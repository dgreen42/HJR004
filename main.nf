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
	"""
}

process index_ref {
	label "HJR004"
	cpus 4
	memory "16 GB"
	input:
		path ref
	output:
		path "index.mmi", emit: index
	script:
	"""
		minimap2 -t 4 -d "index.mmi" ${ref}
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
		path "${fastq.baseName}_aln.bam", emit: bampath
	script:
	"""
		minimap2 -ax splice -uf -k14 ${idx} ${fastq} \
		| samtools view -q 10 -m 20 -F 0x304 -o ${fastq.baseName}_aln.bam
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
		val bampath
	output:
		tuple file("bambu_out/HJR004_counts_gene.txt"), file("bambu_out/HJR004_counts_trasncript.txt"), file("bambu_out/HJR004_CPM_transcript.txt"), file("bambu_out/HJR004_extended_annotations.gtf"), file("bambu_out/HJR004_fullLengthCounts_transcripts.txt"), file("bambu_out/HJR_uniqueCounts_transcript.txt")
	script:
	"""
		bambu_cts.R ${genome} ${anno} ${ndr} ${bampath}
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

	fq_ch = channel.fromPath(params.fastq + "*.fastq.gz")
	check_tools()
	index_ref(ref_anno)
	align_and_map(fq_ch, index_ref.out.index, ref_genome)
	count_transcripts(ref_genome, ref_anno, params.ndr, align_and_map.out.bampath)

}

workflow.onComplete {
	println("All done!")
}
workflow.onError {
	println("Well shit..." + error)
}


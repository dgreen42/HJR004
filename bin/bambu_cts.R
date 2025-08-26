#!/usr/bin/env Rscript
suppressMessages(library(bambu))

args <- commandArgs(trailingOnly = T)

#arg1: genome
#arg2: annotation
#arg3: ndr

print(args[1])
print(args[2])
print(args[3])

genome <- file.path(args[1])
annotation <- file.path(args[2])
prepAnno <- prepareAnnotations(annotation)
samples <- list.files(path = "aln", recursive = T, full.names = T)

analysis <- bambu(reads = samples, annotations = prepAnno, genome = genome, NDR = as.numeric(args[3]))
writeBambuOutput(analysis, path = "bambu_out", prefix = "HJR004_")
#writeToGTF(rowRanges(analysis), file = "./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")


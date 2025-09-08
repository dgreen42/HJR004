#!/usr/bin/env Rscript
library(bambu)
library(tools)

source("../../../bin/utils.R")
args <- commandArgs(trailingOnly = T)

#arg1: genome
#arg2: annotation
#arg3: ndr
#arg4..n: bams

hold_args <- c()
for(i in 4:length(args)) {
	hold_args[i] <- file_path_as_absolute(args[i])
}
samples <- ramna(hold_args)
print(samples)


genome <- file.path(args[1])
annotation <- file.path(args[2])
prepAnno <- prepareAnnotations(annotation)
samples <- hold_args[4:length(args)]
rm(hold_args)

# warnings are suppressed because bambu needs to fix its use of depricated function 'summarise' in dplyr
suppressWarnings(analysis <- bambu(reads = samples, annotations = prepAnno, genome = genome, NDR = as.numeric(args[3])))
writeBambuOutput(analysis, path = "bambu_out", prefix = "HJR004_")


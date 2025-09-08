#!/usr/bin/env Rscript
library(stageR)
library(edgeR)
library(DEXSeq)

source("../../../bin/utils.R")
args <- commandArgs(trailingOnly = T)

#arg1: counts file
#arg2: anno
#arg3: sample_sheet

#file_path <- file_path_as_absolute(args[1])
file_path <- args[1]
sample_sheet_path <- args[3]
print(file_path)
print(sample_sheet_path)
sample_sheet <- read.csv(sample_sheet_path)
cts <- read.delim(file_path)
rownames(cts) <- cts$TXNAME
cts <- cts[,3:ncol(cts)]
sample_sheet
cols <- colnames(cts)
cols
group <- c()
for (i in 1:length(cols)) {
    colsp <- strsplit(cols[i], "[.]")[[1]][1]
    for (j in 1:nrow(sample_sheet)) {
        if (colsp == sample_sheet$sample[j]) {
            group[i] <- trimws(sample_sheet$treatment[j])
        }
    }
}
group
dge <- DGEList(counts = cts, group = group)
head(dge)
expDesign <- model.matrix(~0+group, data = dge$samples)
colnames(expDesign) <- levels(dge$samples$group)
cmpOffset <- 2
keep <- rowSums(cpm(round(cts, digit = 0))>cmpOffset) >= 2
dge <- DGEList(cts[keep,])
colnames(dge) <- rownames(expDesign)
dge <- calcNormFactors(dge)
voomObj <- voom(dge, expDesign, plot=T)
fit <- lmFit(voomObj, expDesign)
contrast.matrix <- makeContrasts(nod-irt, nod-mrt, irt-mrt, levels = expDesign)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
de <- decideTests(fit2)
summary.TestResults(de)

write.csv(fit2$coefficients, file = "de_coefficients.csv")
write.csv(de, file = "de_results.csv")
write.csv(summary.TestResults(de), file = "de_summary.csv")

# stage wise analysis

alpha <- 0.05
nGenes <- nrow(dge)
tableF <- topTableF(fit2, number = nGenes, sort.by = "none")
pScreen <- tableF$P.Value
names(pScreen) = rownames(tableF)
pConf <- sapply(1:3, function(i) topTable(fit2, coef = i, number = nGenes, sort.by = "none")$P.Value)
dimnames(pConf) <- list(rownames(fit2), c("n-i", "n-m", "i-m"))
stageRObj <- stageR(pScreen = pScreen, pConfirmation = pConf, pScreenAdjusted = F)
stageRObj <- stageWiseAdjustment(object = stageRObj, method = "none", alpha = alpha)
res <- getResults(stageRObj)

write.csv(res, "adjusted_results.csv")
write.csv(colSums(res), "adjusted_results_summary.csv")

print("de done")

# differential exon expression ----

cts <- read.delim(args[1])
anno <- read.delim(args[2])

rownames(cts) <- cts$TXNAME

tx2gene <- cts[,1:2]
colnames(tx2gene) <- c("transcript", "gene")
tcts <- cts[,3:ncol(cts)]
tx2dist <- table(table(tx2gene$gene))
write.csv(tx2dist, "taxa_to_gene_distribution.csv")
sampleData <- read.csv(args[3])
rownames(sampleData) <- colnames(tcts)
sampleData <- sampleData[,2:3]
# this is the same as cts$GENEID
geneForEachTx <- tx2gene[match(rownames(tcts), tx2gene[,1]),2]
dxd <- DEXSeqDataSet(countData = round(tcts, digit = 0),
                     sampleData = sampleData,
                     design = ~ sample + exon + group + treatment:exon,
                     featureID = rownames(tcts),
                     groupID = as.character(geneForEachTx)
)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon + group)
dxr <- DEXSeqResults(dxd)
qvalDxr <- perGeneQValue(dxr)


# stage wise analysis

pConf <- matrix(dxr$pvalue, ncol = 1)
dimnames(pConf) <- list(c(dxr$featureID), c("transcript"))
pScreen <- qvalDxr
genesub <- rep(NA, length(tx2gene$gene))
n <- 0
for(gene in tx2gene$gene) {
    n <- n + 1
    genesub[n] <- gsub(" ", "", gene)
}
tx2gene$gene <- genesub
stageRTxObj <- stageRTx(pScreen = pScreen,
                        pConfirmation = pConf,
                        pScreenAdjusted = T,
                        tx2gene = tx2gene
                        )
stageRTxObj <- stageWiseAdjustment(object = stageRTxObj, method = "dtu", alpha = 0.05, allowNA = T)
padj <- getAdjustedPValues(stageRTxObj, order = T, onlySignificantGenes = T)
props <- isoformProp2(cts)

altcolnames <- c("start", "end", "transcripts", "acronym", "strand")
altsplice <- data.frame(matrix(NA, nrow = 1, ncol = 5))
colnames(altsplice) <- altcolnames
count <- 0

for(gene in unique(padj$geneID)) {
    count <- count + 1
    altsplice[count,] <- plotIsoform(strsplit(gene, ";")[[1]][2], annotation = "./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_extended_annotations.gtf" , exon_marker = T,
            acronym_list = "./LeGOO-Currated-GeneNames.csv",
	    suppress_plot = T)
}

write.csv(padj, "dex_adjusted_pval.csv")
write.csv(props, "dex_isoform_proportions.csv")
write.csv(altsplice, "dex_altsplice.csv")

library(stageR)
library(edgeR)

cts <- read.delim("./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_counts_transcript.txt")
rownames(cts) <- cts$TXNAME
cts <- cts[,3:ncol(cts)]
head(cts)
group <- rep(c("nod" ,"irt" ,"mrt"), 3)
dge <- DGEList(counts = cts, group = group)
head(dge)
expDesign <- model.matrix(~0+group, data = dge$samples)
colnames(expDesign) <- levels(dge$samples$group)
expDesign
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
colSums(res)

resCounts <- count_udns(res)

write.csv(, file = "./bambu_out_NDR_3_filtered/stageR_significant_splice/differential_expr.csv", row.names = F)

# differential exon expression ----
library(DEXSeq)

cts <- read.delim("./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_counts_transcript.txt")
anno <- read.delim("./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", sep = "\t", header = F)

rownames(cts) <- cts$TXNAME

tx2gene <- cts[,1:2]
colnames(tx2gene) <- c("transcript", "gene")
tcts <- cts[,3:ncol(cts)]
tx2dist <- table(table(tx2gene$gene))
barplot(tx2dist, main = "Dist of Transcript per Gene")
sampleData <- read.csv("./sample_data.csv")
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

write.csv(padj, file = "./bambu_out_NDR_3_filtered/stageR_significant_splice/adjusted_pval.csv", row.names = F)

source("utils.R")

padj <- read.csv("./bambu_out_NDR_3_filtered/stageR_significant_splice/adjusted_pval.csv")

altcolnames <- c("start", "end", "transcripts", "acronym", "strand")
altsplice <- data.frame(matrix(NA, nrow = 1, ncol = 5))
colnames(altsplice) <- altcolnames
count <- 0

for(gene in unique(padj$geneID)) {
    count <- count + 1
    altsplice[count,] <- plotIsoform(strsplit(gene, ";")[[1]][2], annotation = "./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_extended_annotations.gtf" , exon_marker = T,
            acronym_list = "./LeGOO-Currated-GeneNames.csv")
}

SYP132 <- plotIsoform("MtrunA17_Chr2g0322801",
            annotation = "./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_extended_annotations.gtf" ,
            exon_marker = T,
            acronym_list = "./LeGOO-Currated-GeneNames.csv"
            )

plotIsoform("MtrunA17_Chr8g0390331",
            annotation = "./bambu_out_NDR_3_filtered/Harris_JR_RNA_004_NDR_3_extended_annotations.gtf" ,
            exon_marker = T,
            acronym_list = "./LeGOO-Currated-GeneNames.csv"
            )

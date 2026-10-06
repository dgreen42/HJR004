# These are utility functions
library(RColorBrewer)
library(foreach)

tryAcroGrep <- function(gene, acronym_list) {
    tryCatch(acro_grep <- system2("grep", args = paste(gene, acronym_list), stdout = T),
             error = function(cond) {
                 message(paste(gene, "does not have an acronym yet"))
                 message(conditionMessage(cond))
             },
             warning = function(cond) {
                 message(conditionMessage(cond))
             },
             finally = {
                 if (exists("acro_grep")) {
                     return(acro_grep)
                 } else {
                     return(0)
                 }
             }
    )
}

plotIsoform <- function(gene, annotation, exon_marker = F, prop = NULL, acronym_list = NULL, conversion_table = NULL, cts = NULL, suppress_plot = F) {
    par(xpd = F)
    tryCatch(grepd <- system2("grep", args = paste(gene, annotation), stdout = T),
             error = function(cond) {
                 message(conditionMessage(cond))
                 return(list(start = NA, 
                             end = NA, 
                             parent_transcript = NA,
                             transcripts = NA,
                             acronym = NA,
                             strand = NA))
             },
             warning = function(cond) {
                 message(conditionMessage(cond))
             },
             finally = {
                 if (exists("grepd")) {
                     message("grep successful")
                 } else {
                     return(list(start = NA, 
                                 end = NA,
                                 parent_transcript = NA,
                                 transcripts = NA,
                                 acronym = NA,
                                 strand = NA))
                 }
             }
    )
    
    rawFeatures <- strsplit(grepd, split = "\t")
    featureFrame <- data.frame(matrix(NA, ncol = length(rawFeatures[[1]]), nrow = length(rawFeatures)))
    colnames(featureFrame) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    for (i in 1:length(rawFeatures)) {
        featureFrame[i,] <- rawFeatures[[i]]
    }
    featureFrame$attribute <- strsplit(featureFrame$attribute, split = ";")
    
    for (i in 1:nrow(featureFrame)) {
        if (featureFrame[i,]$feature == "transcript") {
            geneid <- paste(strsplit(featureFrame[i,]$attribute[[1]][1], split = "\"")[[1]][2],
                            gsub("\"", "", featureFrame[i,]$attribute[[1]][2]))
            transcriptid <- strsplit(featureFrame[i,]$attribute[[1]][3], split = "\"")[[1]][2]
            featureFrame$geneid[i] = geneid
            featureFrame$transcriptid[i] = transcriptid 
            featureFrame$exonid[i] = NA
        } else if (featureFrame[i,]$feature == "exon") {
            geneid <- paste(strsplit(featureFrame[i,]$attribute[[1]][1], split = "\"")[[1]][2],
                            gsub("\"", "", featureFrame[i,]$attribute[[1]][2]))
            transcriptid <- strsplit(featureFrame[i,]$attribute[[1]][3], split = "\"")[[1]][2]
            exonid <- strsplit(featureFrame[i,]$attribute[[1]][4], split = "\"")[[1]][2]
            featureFrame$geneid[i] = geneid
            featureFrame$transcriptid[i] = transcriptid 
            featureFrame$exonid[i] = exonid 
        }
    }
    
    gene_name <- NULL
    if (is.null(acronym_list)) {
        gene_name <- gene
    } else {
        acro_grep <- tryAcroGrep(gene, acronym_list = acronym_list)
        if (acro_grep != 0) {
            gene_name <- strsplit(acro_grep, ",")[[1]][2]
        } else {
            gene_name <- gene
        }
    }
    
    if (suppress_plot == T) {
        prop <- NULL
    } else {
        prop <- prop
    }
    
    if (is.null(prop) & suppress_plot != T) {
        transcripts <- unique(featureFrame$transcriptid)
        if(!is.null(conversion_table)) {
            transcript_names <- convertIsoformNames(conversion_table, transcripts,
                                               paste("gene_biotype mRNA;", gene)
                                               )
            if(is.na(transcript_names[1])) {
                transcript_names <- transcripts
            }
        }
        par(mai = c(1.02,2,1,0.42))
        xlimit =  c(as.integer(min(featureFrame$start)), as.integer(max(featureFrame$end)))
        if (is.null(gene_name)) {
            plot(NA,
                 xlim = xlimit,
                 ylim = c(0,length(transcripts)),
                 main = gene,
                 xlab = "Location (bp)",
                 ylab = NA,
                 yaxt = "n",
                 bty = "n",
            )
        } else {
            plot(NA,
                 xlim = xlimit,
                 ylim = c(0,length(transcripts)),
                 main = gene_name,
                 xlab = "Location (bp)",
                 ylab = NA,
                 yaxt = "n",
                 bty = "n",
            )
        } 
        axis(side = 2,
             at = 1:length(transcripts),
             labels = transcript_names,
             las = 1,
             cex.axis = 0.8)
        count <- 0
        for (i in transcripts) {
            count <- count + 1
            transFrame <- featureFrame[featureFrame$transcriptid == i, ]
            for(j in 1:nrow(transFrame)) {
                row <- transFrame[j,]
                if (row$feature == "transcript") {
                    lines(x = c(row$start, row$end), y = c(count,count), lwd = 1, lend = 1)
                } else if (row$feature == "exon") {
                    lines(x = c(row$start, row$end), y = c(count,count), lwd = 10, lend = 1)
                    if (exon_marker == T) {
                        abline(v = row$start, lty = 2)
                        abline(v = row$end, lty = 2)
                    }
                }
            }
        }
        mtext("Isoform",
              side = 3,
              padj = -1.2, 
              adj = -0.3,
        )
        mtext(paste("Strand:", featureFrame$strand[1]),
              side = 3,
              padj = -1.2,
              adj = -0.01
        )
        ts <- ""
        for(i in transcripts) {
            if (i != gene) {
                ts <- paste(ts, i)
            }
        }
        return(list(start = xlimit[1], 
                    end = xlimit[2], 
                    parent_transcript = gene,
                    alternative_transcripts = trimws(ts),
                    acronym = gene_name,
                    strand = featureFrame$strand[1])) 
    } else if (suppress_plot != T) {
        transcripts <- unique(featureFrame$transcriptid)
        if(!is.null(conversion_table)) {
            transcript_names <- convertIsoformNames(conversion_table, transcripts,
                                               paste("gene_biotype mRNA;", gene)
                                               )
            if(is.na(transcript_names[1])) {
                transcript_names <- transcripts
            }
        }
        par(mai = c(1.02,2,1,1))
        xlimit =  c(as.integer(min(featureFrame$start)), as.integer(max(featureFrame$end)))
        if (is.null(gene_name)) {
            plot(NA,
                 xlim = xlimit,
                 ylim = c(0,length(transcripts)),
                 main = gene_name,
                 xlab = "Location (bp)",
                 ylab = NA,
                 yaxt = "n",
                 bty = "n",
            )
        } else {
            plot(NA,
                 xlim = xlimit,
                 ylim = c(0,length(transcripts)),
                 main = gene_name,
                 xlab = "Location (bp)",
                 ylab = NA,
                 yaxt = "n",
                 bty = "n",
            )
        }
        axis(side = 2,
             at = 1:length(transcripts),
             labels = transcript_names,
             las = 1,
             cex.axis = 0.8)
        count <- 0
        for (i in transcripts) {
            count <- count + 1
            transFrame <- featureFrame[featureFrame$transcriptid == i, ]
            for(j in 1:nrow(transFrame)) {
                row <- transFrame[j,]
                if (row$feature == "transcript") {
                    lines(x = c(row$start, row$end), y = c(count,count), lwd = 1, lend = 1)
                } else if (row$feature == "exon") {
                    lines(x = c(row$start, row$end), y = c(count,count), lwd = 10, lend = 1)
                    if (exon_marker == T) {
                        abline(v = row$start, lty = 2)
                        abline(v = row$end, lty = 2)
                    }
                }
            }
        }
        propidx <- grep(gene, prop$GENEID)
        props <- c()
        plist <- data.frame(matrix(NA, ncol = ncol(prop)))
        colnames(plist) <- colnames(prop)
        count <- 1
        for(i in propidx) {
            plist[count,] <- prop[i,]
            count <- count + 1
        }
        count <- 1
        for (i in transcripts) {
            for (j in 1:nrow(plist)) {
                if (i == plist$TXNAME[j]) {
                    props[count] <- round(as.double(plist$prop[j]), digit = 4)
                    count <- count + 1
                } else {
                    next
                }
            }
        }
        if (length(props) > 4) {
            if (length(props) >= 6) {
                font.size = 0.8
            } else {
                font.size = 1
            }
        } else {
            font.size = 1
        }
        axis(side = 4,
             at = 1:length(propidx),
             labels = props,
             cex.axis = font.size,
             las = 1
        )
        mtext("Isoform Proportion",
              side = 3,
              adj = 1.2,
              padj = -1.2
        )
        mtext("Isoform",
              side = 3,
              padj = -1.2, 
              adj = -0.2,
              
        )
        mtext(paste("Strand:", featureFrame$strand[1]),
              side = 3,
              padj = -1.2,
              adj = -0.01
        )
        ts <- ""
        for(i in transcripts) {
            if (i != gene) {
                ts <- paste(ts, i)
            }
        }
        return(list(start = xlimit[1],
                    end = xlimit[2],
                    parent_transcript = gene,
                    alternative_transcripts = trimws(ts),
                    acronym = gene_name,
                    strand = featureFrame$strand[1],
                    props = plist
                    )
        )
    } else { 
        transcripts <- unique(featureFrame$transcriptid)
        if(!is.null(conversion_table)) {
            transcript_names <- convertIsoformNames(conversion_table, transcripts,
                                               paste("gene_biotype mRNA;", gene)
                                               )
            if(is.na(transcript_names[1])) {
                transcript_names <- transcripts
            }
        }
        ts <- ""
        for(i in transcripts) {
            ts <- paste(ts, i)
        }
        xlimit =  c(as.integer(min(featureFrame$start)), as.integer(max(featureFrame$end)))
        return(list(start = xlimit[1], 
                    end = xlimit[2], 
                    parent_transcript = gene,
                    transcripts = ts,
                    acronym = gene_name,
                    strand = featureFrame$strand[1])) 
    }
}

plotVen <- function(left, center, right, title, labl = NULL, labc = NULL, labr = NULL) {
    par(xpd = NA)
    plot(NA, xlim = c(0,100), ylim = c(0,100),
         xlab = NA, ylab = NA,
         xaxt = "n", yaxt = "n",
         bty = "n")
    points(35, 40, cex = 40, col = rgb(252, 211, 3, 150, maxColorValue = 255), pch = 16)
    points(65, 40, cex = 40, col = rgb(3, 78, 252, 150, maxColorValue = 255), pch = 16)
    text(20, 40, paste(left), cex = 2)
    text(50, 40, paste(center), cex = 2)
    text(80, 40, paste(right), cex = 2)
    title(paste(title))
    text(20, -30, paste(labl), cex = 1.5)
    text(50, -30, paste(labc), cex = 1.5)
    text(80, -30, paste(labr), cex = 1.5)
}

ramna <- function(x) {
    y <- NULL
    count <- 1
    for (i in x) {
        if (is.na(i)) {
            next
        } else {
            y[count] <- i
            count <- count + 1
        }
    }
    return(y)
}

nan.to.zero <- function(x) {
    y <- NULL
    count <- 1
    for (i in x) {
        if (is.nan(i)) {
            y[count] <- 0
            count <- count + 1
        } else {
            y[count] <- i
            count <- count + 1
        }
    }
    return(y)

}

compareReg <- function(set1, set2) {
    summary <- list(
        upReg = NA,
        downReg = NA,
        noSig = NA
    )
    countup <- 0
    countdown <- 0
    countnosig <- 0
    reg <- list(
        up = c(),
        down = c(),
        summary = list()
    )
    for (i in set1$upReg) {
        for(j in set2$upReg) {
            if (i == j) {
                countup <- countup + 1
                reg$up[countup] <- i
            }
        }
    }
    
    summary$upReg <- countup
    
    for (i in set1$downReg) {
        for(j in set2$downReg) {
            if (i == j) {
                countdown <- countdown + 1
                reg$down[countdown] <- i
            }
        }
    }
    summary$downReg <- countdown
    
    for (i in set1$noSig) {
        for(j in set2$noSig) {
            if (i == j) {
                countnosig <- countnosig + 1
            }
        }
    }
    summary$noSig <- countnosig
    reg$summary <- summary
    
    return(reg)
}

getUniqueReg <- function(reg, compReg) {
    regUpUnique <- c()
    count <- 1
    for (i in reg) {
        test <- sum(compReg == i)
        if (test == 1) {
            next
        } else {
            regUpUnique[count] <- i
            count <- count + 1
        }
    }
    return(regUpUnique)
}

getGeneIDsRegUnique <- function(regUnique, master, geneLookup, n) {
    if (length(regUnique) != 0) {
        count <- 1
        for(i in regUnique) {
            for(j in 1:nrow(geneLookup)) {
                if (i == geneLookup$TXNAME[j]) {
                    master[[n]][count,] <- geneLookup[j,]
                    count <- count + 1
                } else {
                    next
                }
            }
        }
    } else {
        print("Comp Reg has no length")
    }
    return(master)
}

regUnique <- function(reg1, reg2, compReg, master, geneLookup, setNames = NULL) {
    print("Starting Reg1")
    regUpUnique1 <- getUniqueReg(reg1, compReg)
    print("Reg1 complete")
    print("Starting Reg2")
    regUpUnique2 <- getUniqueReg(reg2, compReg)
    print("Reg2 complete")
    
    print("Starting Gene Lookup 1")
    master <- getGeneIDsRegUnique(regUpUnique1, master, geneLookup, 1)
    print("Lookup 1 complete")
    print("Starting Lookup 2")
    master <- getGeneIDsRegUnique(regUpUnique2, master, geneLookup, 2)
    print("Lookup 2 complete")
    print("Starting Comp Lookup")
    master <- getGeneIDsRegUnique(compReg, master, geneLookup, 3)
    print("Comp Lookup Complete")
    
    if (!is.null(setNames)) {
        names(master) <- setNames
    }
    
    return(master)
}

isoformProp2 <- function(counts) {
    t <- 0
    f <- 0
    s <- 0
    len <- nrow(counts)

    props <- data.frame(TXNAME = NA, GENEID = NA, genetotal = NA, isototal = NA, prop = NA)
    nodprop <- data.frame(TXNAME = NA, GENEID = NA, genetotal = NA, isototal = NA, prop = NA)
    irtprop <- data.frame(TXNAME = NA, GENEID = NA, genetotal = NA, isototal = NA, prop = NA)
    mrtprop <- data.frame(TXNAME = NA, GENEID = NA, genetotal = NA, isototal = NA, prop = NA)

    genes <- unique(counts$GENEID)
    count <- 1

    nod <- sampleData[sampleData$group == "nod"]
    irt <- sampleData[sampleData$group == "irt"]
    mrt <- sampleData[sampleData$group == "mrt"]

    nodidx <- getColIdx(cts, nod)
    irtidx <- getColIdx(cts, irt)
    mrtidx <- getColIdx(cts, mrt)

    nodsub <- createSubsetCts(cts, nodidx)
    irtsub <- createSubsetCts(cts, irtidx)
    mrtsub <- createSubsetCts(cts, mrtidx)

    head(nodsub)
    head(irtsub)
    head(mrtsub)

    for(gene in genes) {
        set <- counts[counts$GENEID == gene,]
        nodset <- nodsub[nodsub$GENEID == gene,]
        irtset <- irtsub[irtsub$GENEID == gene,]
        mrtset <- mrtsub[mrtsub$GENEID == gene,]
        genetotal <- sum(set[,3:ncol(set)])
        nodgenetotal <- sum(nodset[,3:ncol(nodset)])
        irtgenetotal <- sum(irtset[,3:ncol(irtset)])
        mrtgenetotal <- sum(mrtset[,3:ncol(mrtset)])

        for(i in 1:nrow(set)) {
            props[count,] <- c(set$TXNAME[i],
                               set$GENEID[i],
                               genetotal,
                               sum(set[i,3:ncol(set)]),
                               sum(set[i,3:ncol(set)])/genetotal
            )

            nodprop[count,] <- c(nodset$TXNAME[i],
                               nodset$GENEID[i],
                               nodgenetotal,
                               sum(nodset[i,3:ncol(nodset)]),
                               sum(nodset[i,3:ncol(nodset)])/nodgenetotal
            )

            irtprop[count,] <- c(irtset$TXNAME[i],
                               irtset$GENEID[i],
                               irtgenetotal,
                               sum(irtset[i,3:ncol(irtset)]),
                               sum(irtset[i,3:ncol(irtset)])/irtgenetotal
            )

            mrtprop[count,] <- c(mrtset$TXNAME[i],
                               mrtset$GENEID[i],
                               mrtgenetotal,
                               sum(mrtset[i,3:ncol(mrtset)]),
                               sum(mrtset[i,3:ncol(mrtset)])/mrtgenetotal
            )

            count <- count + 1
            if (count > len*0.25 && count < len*0.251 && t == 0) {
                print("--- 25% complete ---")
	    	t <- 1
            } else if (count > len*0.5 && count < len*0.51 && f == 0) {
                print("--- 50% complete ---")
	    	f <- 1
            } else if (count > len*0.75 && count < len*0.751 && s == 0) {
                print("--- 75% complete ---")
	    	s <- 1
            }
        }
    }

    write.csv(props, "dex_isoform_proportions.csv")
    write.csv(nodprop, "dex_isoform_proportions_nod.csv")
    write.csv(irtprop, "dex_isoform_proportions_irt.csv")
    write.csv(mrtprop, "dex_isoform_proportions_mrt.csv")
}

isoformProp3 <- function(counts, padj) {
    cores <- parallel::detectCores() - 1
    cluster <- parallel::makeCluster(
        cores,
	type = "PSOCK"
    )
    doParallel::registerDoParallel(cl = cluster)
    print(foreach::getDoParRegistered())
    print(foreach::getDoParWorkers())

    genes <- unique(padj$geneID)
    isoprops <- foreach(
        gene = genes,
        .combine = "rbind"
    ) %do% {
        count <- 1
        set <- counts[gsub(" ", "", counts$GENEID) == gene,]
        genetotal <- sum(set[,3:ncol(set)])
        prop <- data.frame(TXNAME = NA, GENEID = NA, genetotal = NA, isototal = NA, prop = NA)
        for(i in 1:nrow(set)) {
            prop[count,] <- c(set$TXNAME[i],
                              set$GENEID[i],
                              genetotal,
                              sum(set[i,3:ncol(set)]),
                              sum(set[i,3:ncol(set)])/genetotal
            )
            count <- count + 1
        }
        prop
    }
    return(isoprops)
}

plotPropComp <- function(propTable, gene) {
    subset <- propTable[propTable$GENEID == gene,]
    complist <- c()
    count <- 1
    for(i in 1:nrow(subset)) {
        for (j in 3:ncol(subset)) {
            complist[count] <-subset[i,j]
            count <- count + 1
        }
    }
    div <- length(complist) / nrow(subset)
    cols <- c()
    for(i in 1:nrow(subset)) {
        cols <- append(cols, rep(colors()[i*4], div))
    }
    barplot(complist,
            width = 1,
            cex.names = 0.6,
            col = cols
    )
    namelocs <- c()
    for(i in 1:nrow(subset)) {
        namelocs[i] <- div * i 
    }
    axis(1,
         at = namelocs,
         labels = subset$TXNAME,
         cex.axis = 0.5
    )
}

plotPropCompSep <- function(propTable, gene) {
    subset <- propTable[propTable$GENEID == gene,]
    layout(matrix(1:4, ncol = 4, nrow = 1))
    par(mai = c(1.3, 0.5, 0.8, 0.5))
    title("Transcript Proportions Across Tissues")
    barplot(subset$NOD, main = "NOD", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
    barplot(subset$IRT, main = "IRT", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
    barplot(subset$MRT, main = "MRT", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
    barplot(subset$all, main = "ALL", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
}

plotPropHeatamp <- function(propTable, gene) {
    subset <- propTable[propTable$GENEID == gene,]
    propMat <- matrix(0, nrow = nrow(subset), ncol = ncol(subset) - 2)
    for (i in 1:nrow(subset)) {
        for (j in 3:ncol(subset)) {
            propMat[i,j - 2] <- as.numeric(subset[i,j])
        }
    }
    rownames(propMat) <- subset$TXNAME
    colnames(propMat) <- colnames(subset[,3:ncol(subset)])
    n <- 10
    nums <- round(seq(0, 1, length.out = n), digit = 1)
    cols <- hcl.colors(n, palette = "Reds")
    hm <- heatmap(
        propMat, Rowv = NA, Colv = NA, col = cols,
        cexRow = 0.8,
        cexCol = 0.8,
    )
    legend("bottomleft",
           legend = nums,
           col = c(cols),
           pch = 15,
           pt.cex = 3
    )
    return(hm)
}

coolMap <- function(data, cluster = F) {
    m <- as.matrix(data)
    d <- dist(m)
    h <- hclust(d)
    par(xpd = T, mai = c(1.3, 1, 0.4, 3.5))
    plot(NA, xlim = c(0, ncol(m)), ylim = c(0, nrow(m)),
         xaxt = 'n',
         yaxt = 'n',
         bty = 'n',
         xlab = NA,
         ylab = NA
    )
    color <- NULL
    
    if(cluster == T) {
        for(i in 1:ncol(m)) {
            count <- 1
            for(j in h$order) {
                if(m[j,i] == 1) {
                    color <- "red"
                } else if (m[j,i] == -1) {
                    color <- "blue"
                } else if (m[j,i] == 0) {
                    color <- "white"
                } else {
                    print("what?")
                }
                polygon(
                    border = NA,
                    x = c(i-1, i, i, i-1),
                    y = c(count-1, count-1, count, count),
                    col = color
                )
                count <- count + 1
            }
            legend("right",
                   legend = c("Upregulated", "No signigicant change", "Downregulated"),
                   col = c("red", "white", "blue"),
                   pch = 15,
                   inset = c(-1, 0)
            )
        }
    } else {
        for(i in 1:ncol(m)) {
            for(j in 1:nrow(m)) {
                if(m[j,i] == 1) {
                    color <- "red"
                } else if (m[j,i] == -1) {
                    color <- "blue"
                } else if (m[j,i] == 0) {
                    color <- "black"
                } else {
                    print("what?")
                }
                polygon(
                    border = NA,
                    x = c(i-1, i, i, i-1),
                    y = c(j-1, j-1, j, j),
                    col = color
                )
                color <- NULL
            }
            legend("right",
                   legend = c("Upregulated", "No signigicant change", "Downregulated"),
                   col = c("red", "black", "blue"),
                   pch = 15,
                   inset = c(0, 0)
            )
        }
    }
    axis(1,
         at = c(0.5,1.5,2.5),
         labels = colnames(data),
         las = 2
    )
}

plotHeatmapIso <- function(propTable, gene, ctable = NULL, color, title = NULL, include_total = F, order = NULL) {
    if (!is.data.frame(propTable)) {
        stop("must provide a data.frame")
    } else {
        subset <- propTable[propTable$GENEID == gene,]
        propMat <- matrix(0, nrow = nrow(subset), ncol = ncol(subset) - 2)
        for (i in 1:nrow(subset)) {
            for (j in 3:ncol(subset)) {
                propMat[i,j - 2] <- as.numeric(subset[i,j])
            }
        }
        
        rownames(propMat) <- subset$TXNAME
        if(!is.null(ctable)) {
            converted_names <- c()
            for(name in rownames(propMat)) {
                cname <- ctable[ctable$BAMBUTX == name,]$acronym
                converted_names <- append(converted_names, cname)
            }
            rownames(propMat) <- converted_names
        }
        colnames(propMat) <- colnames(subset[,3:ncol(subset)])
        
        if(!isTRUE(include_total)) {
            propMat <- propMat[,1:3]
        }
        
        divx <- seq(1, ncol(propMat))
        divy <- seq(1, nrow(propMat))
        xlimit = c(0, ncol(propMat))
        ylimit = c(0, nrow(propMat))
        par(xpd = T, mai = c(0.5, 2, 0.3, 1.5))
        plot(NA, xlim = xlimit, ylim = ylimit,  bty = "n", xaxt = "n", yaxt = "n",
             xlab = NA, ylab = NA)
        xaxisat <- 1:ncol(propMat) - 0.5
        yaxisat <- 1:nrow(propMat) - 0.5
        axis(1, at = xaxisat, labels = colnames(propMat), padj = -2, tick = F)
        axis(2, at = yaxisat, labels = rownames(propMat), hadj = 0.85, tick = F, las = 2)
        propMat100 <- round(propMat * 100, digit = 0) + 1
        colors <- colorRampPalette(color)(max(propMat100))
        for (i in divx) {
            for (j in divy) {
                polygon(
                    x = c(i - 1, i, i, i - 1),
                    y = c(j - 1, j - 1, j, j),
                    col = colors[propMat100[j,i]]
                )
            }
        }
        nums <- round(seq(0, max(propMat), length.out = 10), digit = 2)*100
        lcolnums <- round(seq(1, length(colors), by = length(colors)/length(nums)), digit = 0)
        lcol <- c()
        count <- 1
        for(i in lcolnums) {
            lcol[count] <- colors[i]
            count <- count + 1
        }
        xmin <- max(xlimit) + 0.25
        xmax <- max(xlimit) + 1.25
        ymin <- max(ylimit) * 0.2
        ymax <- max(ylimit) * 0.92
        xleg <- c(xmin, xmax, xmax, xmin)
        yleg <- c(ymin, ymin, ymax, ymax)
        polygon(x = xleg,
                y = yleg
        )
        xadd <- 0.25
        yadd <- 0.1
        xsizer <- 0.2
        ysizer <- max(ylimit) *0.0667
        for(i in 1:length(nums)) {
            xcoord <- xmin + xadd
            ycoord <- ymin + yadd
            polygon(x = c(xcoord, xcoord + xsizer, xcoord + xsizer, xcoord),
                    y = c(ycoord, ycoord, ycoord + ysizer, ycoord + ysizer),
                    col = lcol[i]
            )
            text(x = xcoord + xsizer + 0.25,
                 y = ycoord + ysizer - 0.1,
                 labels = nums[i]
            )
            yadd <- yadd + ysizer
        }
        axis(4, 4, "% abundance of sequenced transcripts", tick = F, padj = 6, hadj = 0.8)
        if(!is.null(title)) {
            title(title)
        }
    }
    
    returnl <- list(
        propMat = propMat,
        xlim = xlimit,
        ylim = ylimit,
        xleglim = xleg,
        yleglim = yleg
    )
    
    return(returnl)
}

regTx <- function(DEX) {
    # upReg <- sum(DEX == 1)
    # downReg <- sum(DEX == -1)
    # noSig <- sum(DEX == 0)
    
    upRegTx <- c()
    downRegTx <- c()
    noSigTx <- c()
    for (i in 1:nrow(DEX)) {
        if (DEX[[i]] == 1) {
            upRegTx[i] <- rownames(DEX)[i]
        } else if (DEX[[i]] == -1) {
            downRegTx[i] <- rownames(DEX)[i]
        } else if (DEX[[i]] == 0) {
            noSigTx[i] <- rownames(DEX)[i]
        }
    }
    
    list(
        upReg = ramna(upRegTx),
        downReg = ramna(downRegTx),
        noSig = ramna(noSigTx)
    )
}

initDiagnosicTalbe <- function() {
    return(
        data.frame(
            sample = NULL,
            passRation = NULL,
            meanLen = NULL,
            medianLen = NULL,
            maxLen = NULL,
            minLen = NULL,
            medianQScore = NULL,
            N50 = NULL
        )
    )
}
            
rnaDiagnosicSummary <- function(summary, file) {
    passfilter <- summary$passes_filtering == T
    name <- strsplit(strsplit(file, "/")[[1]][5], "_")[[1]][1]
    summary <- data.frame(
        sample = name,
        passRation = sum(passfilter) / nrow(summary),
        meanLen = mean(summary$sequence_length_template[passfilter]),
        medianLen = median(summary$sequence_length_template[passfilter]),
        maxLen = max(summary$sequence_length_template[passfilter]),
        minLen = min(summary$sequence_length_template[passfilter]), 
        medianQScore = round(median(summary$mean_qscore_template[passfilter]), digits = 2),
        N50 = getNstat(summary, passfilter, 50)
        )
    return(summary)
}

getNstat <- function(summary, filter, N) {
    lens <- as.numeric(sort(summary$sequence_length_template[passfilter], decreasing = T))
    idx <- which(cumsum(lens) / sum(lens) >= 50/100)[1]
    lens[idx]
}

#udns:
#   up regulate
#   down regulated
#   no siginificance
count_udns <- function(results) {
    resultsFrame <- data.frame(matrix(0, nrow = 3, ncol = ncol(results) - 1))
    colnames(resultsFrame) <- colnames(results)[2:ncol(results)]
    rownames(resultsFrame) <- c("Down", "NoSig", "Up")
    for(row in 1:nrow(results)) {
        if (row[1] == 1) {
            for(i in 2:ncol(results)) {
                if (row[i] == 1) {
                    resultsFrame[3,i] <- resultsFrame[3,i] + 1
                } else if (row[i] == -1) {
                    resultsFrame[1,i] <- resultsFrame[1,i] + 1
                } else if (row[i] == 0) {
                    resultsFrame[2,i] <- resultsFrame[2,i] + 1
                }
            }
        } else {
            next
        }
    }
    return(resultsFrame)
}

getColIdx <- function(cts, group) {
    idx <- NULL
    count <- 0
    countidx <- 0
    for(i in colnames(cts)) {
        count <- count + 1
        for(j in rownames(group)) {
            if (i == j) {
                countidx <- countidx + 1
                idx[countidx] <- count
            }
        }
    }
    names(idx) <- rownames(group)
    return(idx)
}

createSubsetCts <- function(cts, subset) {
    subdf <- data.frame(matrix(NA, ncol = length(subset) + 2, nrow = nrow(cts)))
    cnames <-  c("TXNAME", "GENEID", names(subset))
    colnames(subdf) <- cnames
    subdf$TXNAME <- cts$TXNAME
    subdf$GENEID <- cts$GENEID
    #start at 2 since the first two coluns will be txname and geneid
    count <- 2
    for(i in subset) {
        count <- count + 1
        subdf[,count] <- cts[,i]
    }
    
    return(subdf)
}

propsAsNumeric <- function(props) {
    for(i in 3:ncol(props)) {
        props[,i] <- as.numeric(props[,i])
    }
    return(props)
}

copyDFscructure <- function(df, all = F) {
    if (all == T) {
        structure <- data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
    } else if(all == F) {
        structure <- data.frame(matrix(NA, nrow = 1, ncol = ncol(df)))
    } else {
        print("Argument for 'all' must be TRUE or FALSE")
    }
    colnames(structure) <- colnames(df)
    return(structure)
}

filterPScreen <- function(de, screen) {
    screened <- copyDFscructure(de, all = F)
    colnames(screened) <- colnames(de)
    count <- 1
    for(i in 1:nrow(screen)) {
        if(screen$padjScreen[i] != 0) {
            screened[count,] <- de[i,]
            rownames(screened)[count] <- rownames(de)[i]
            count <- count + 1
        }
    }
    return(screened)
}

plotTable <- function(dataframe) {
    if (!is.data.frame(dataframe)) {
        stop("must provide a data.frame")
    } else {
        df <- matrix(0, nrow = nrow(dataframe), ncol = ncol(dataframe))
        for (i in 1:nrow(dataframe)) {
            for (j in 1:ncol(dataframe)) {
                df[i,j] <- as.numeric(dataframe[i,j])
            }
        }
        rownames(df) <- rownames(dataframe)
        colnames(df) <- colnames(dataframe)
        divx <- seq(1, ncol(df))
        divy <- seq(1, nrow(df))
        xlimit = c(0, ncol(df))
        ylimit = c(0, nrow(df))
        par(xpd = T, mai = c(0.5, 2, 1, 1))
        plot(NA, xlim = xlimit, ylim = ylimit,  bty = "n", xaxt = "n", yaxt = "n",
             xlab = NA, ylab = NA)
        xaxisat <- 1:ncol(df) - 0.5
        yaxisat <- 1:nrow(df) - 0.5
        axis(1, at = xaxisat, labels = colnames(df), padj = -2, tick = F)
        axis(2, at = yaxisat, labels = rownames(df), hadj = 0.85, tick = F, las = 2)
        for(i in divx) {
            for(j in divy) {
                polygon(
                    x = c(i - 1, i, i, i - 1),
                    y = c(j - 1, j - 1, j, j),
                    text(i - 0.5, j - 0.5, df[j,i])
                )
            }
        }
    }
}

getDEacro <- function(degs, counts, acrotable) {
    t <- copyDFscructure(acrotable, all = F)
    count <- 1
    for(i in 1:nrow(degs)) {
        gene <- trimws(strsplit(counts$GENEID[counts$GENEID == rownames(degs)[i]], ";")[[1]][2])
        if(is.na(gene)) {
            gene <- gcts$GENEID[gcts$GENEID == rownames(degs)[i]]
        } else {
            info <- acrotable[acrotable[,1] == gene,]
            if(nrow(info) != 0) {
                for(i in 1:nrow(info)) {
                    t[count,] <- info[i,]
                    count <- count + 1
                }
            }
        }
    }
    return(t)
}


getDTUacro <- function(dets, counts, acrotable) {
    t <- copyDFscructure(acrotable, all = F)
    t$transcript_name <- NA
    count <- 1
    
    for(i in 1:nrow(dets)) {
        geneid <- counts$GENEID[counts$TXNAME == rownames(dets)[i]]
        gene <- trimws(strsplit(geneid, ";")[[1]][2])
        if(is.na(gene)) {
            next
        } else {
            info <- acrotable[acrotable[,1] == gene,]
            if(nrow(info) != 0) {
                for(j in 1:nrow(info)) {
                    info$transcript_name <- rownames(dets)[i]
                    t[count,] <- info[j,]
                    count <- count + 1
                }
            }
        }
    }
    return(t)
}

findGenes <- function(table, pattern) {
    df <- copyDFscructure(table, all = F)
    count <- 1
    for(i in 1:nrow(table)) {
        grp <- grepl(pattern, table$ACRONYM[i], ignore.case = T)
        if(grp) {
            df[count,] <- table[i,]
            count <- count + 1
        } else {
            next
        }
    }
    return(df)
}

findGenes2 <- function(table, pattern) {
    df <- copyDFscructure(table, all = F)
    rnames <- NULL
    count <- 1
    for(i in 1:nrow(table)) {
        grp <- grepl(pattern, rownames(table)[i], ignore.case = T)
        if(grp) {
            df[count,] <- table[i,]
            rnames[count] <- rownames(table)[i]
            count <- count + 1
        } else {
            next
        }
    }
    rownames(df) <- rnames
    return(df)
}

plotIsoformLFC <- function(props1, props2, contrast1, contrast2, lfc, gene, acrotable) {
    subset1 <- props1[props1$GENEID == gene,]
    subset2 <- props2[props2$GENEID == gene,]
    
    lfcsubset1 <- copyDFscructure(lfc, all = F)
    lfcsubset2 <- copyDFscructure(lfc, all = F)
    
    propsubset1 <- copyDFscructure(subset1, all = F)
    propsubset2 <- copyDFscructure(subset2, all = F)
    
    count <- 1
    
    for(i in 1:nrow(subset1)) {
        for(j in 1:nrow(lfc)) {
            if(subset1$TXNAME[i] == lfc$TXNAME[j]) {
                lfcsubset1[count,] <- lfc[j,]
                propsubset1[count,] <- subset1[i,]
                
                lfcsubset2[count,] <- lfc[j,]
                propsubset2[count,] <- subset2[i,]
                
                count <- count + 1
            }
        }
    }
    
    lfcsubcont1 <- data.frame(TXNAME = lfcsubset1$TXNAME, CONTRAST = eval(parse(text = paste("lfcsubset1$", contrast1, sep = ""))))
    lfcsubcont2 <- data.frame(TXNAME = lfcsubset2$TXNAME, CONTRAST = eval(parse(text = paste("lfcsubset2$", contrast2, sep = ""))))
    
    par(mai = c(2,0.5,1,2.5), font = 2, xpd = F)
    plot(NA, xlim = c(0,1), ylim = c(1, nrow(lfcsubcont1) + 2),
         xlab = "Isoform Proportion", ylab = NA, yaxt = "n", bty = "n")
    for(i in 1:nrow(lfcsubcont1)) {
        x <- propsubset1$prop[i]
        x2 <- propsubset2$prop[i]
        y <- i + 1
        cont1 <- lfcsubcont1$CONTRAST[i]
        cont2 <- lfcsubcont2$CONTRAST[i]
        lines(c(x, x), c(y, 1), lty = 2)
        lines(c(x, 1), c(y, y), lty = 2)
        lines(c(x2, x2), c(y, 1), lty = 2)
        lines(c(x2, 1), c(y, y), lty = 2)
        if(cont1 > 0) {
            lines(c(x, x), c(y, y + cont1 / 2) , lwd = 3, col = "red")
        } else {
            lines(c(x, x), c(y, y + cont1 / 2) , lwd = 3, col  = "blue")
        }
        if(cont2 > 0) {
            lines(c(x2, x2), c(y, y + cont2 / 2) , lwd = 3, col = "red")
        } else {
            lines(c(x2, x2), c(y, y + cont2 / 2) , lwd = 3, col  = "blue")
        }
        
        points(x, y, pch = 21, bg = "black")
        points(x2, y, pch = 21, bg = "grey")
        
    }
    axis(4, 1:nrow(lfcsubcont1) + 1,
         labels = lfcsubcont1$TXNAME,
         las = 1
    )
    mtext("Isoform",
          side = 3,
          padj = 0, 
          adj = 1.25,
          cex = 1.2
    )
    locus <- trimws(strsplit(gene, ";")[[1]][2])
    acronym <- acrotable[acrotable$Mt5.0..r1.7..locus.tag == locus,]$ACRONYM
    if(acronym != 0) {
        title(paste("Isoform Proportion and Log2-Fold Change of\n", acronym))
    } else {
        title(paste("Isoform Proportion and Log2-Fold Change of\n", locus))
    }
    name1 <- strsplit(contrast1, "vs")[[1]][1]
    name2 <- strsplit(contrast2, "vs")[[1]][1]
    
    par(xpd = T)
    legend("bottom",
           inset = c(-0.5, -0.9),
           legend = c(paste(name1, "proportion in contrast:", contrast1),
                      paste(name2, "proportion in contrast:", contrast2),
                      "Positive LFC",
                      "Negative LFC"
           ),
           col = c("black", "black", "red", "blue"),
           pt.bg = c("black", "grey"),
           pch = c(21, 21, 15, 15),
           text.width = 0.9
    )
}

plotMultiReg <- function(genes, cts, lfc, acrotable, conversion_table) {
    par(mai = c(1,1,1,3))
    for(gene in genes) {
        subset <- cts[cts$GENEID == gene,]
        if(nrow(subset) == 0) {
            next
        } else {
            sub_lfc <- data.frame(matrix(NA, ncol = ncol(lfc)))
            colnames(sub_lfc) <- colnames(lfc)
            count <- 1
            for(i in 1:nrow(subset)) {
                for(j in 1:nrow(lfc)) {
                    if(subset$TXNAME[i] == lfc$TXNAME[j]) {
                        sub_lfc[count,] <- lfc[j,]
                        count <- count + 1
                    }
                }
            }
            
            rnames <- convertIsoformNames(conversion_table, sub_lfc[,1], gene)
            rownames(sub_lfc) <- rnames
            
            sub_lfc <- sub_lfc[,2:ncol(sub_lfc)]
            sub <- as.matrix(sub_lfc)
            sub_lfc <- matrix(NA, nrow = nrow(sub), ncol = ncol(sub))
            colnames(sub_lfc) <- colnames(sub)
            rownames(sub_lfc) <- rownames(sub)
            for(i in 1:nrow(sub)) {
                for(j in 1:ncol(sub)) {
                    sub_lfc[i,j] <- as.numeric(sub[i,j])
                }
            }
            
            print(sub_lfc)
            
            locus <- trimws(strsplit(gene, ";")[[1]][2])
            acronym <- acrotable[acrotable$Mt5.0..r1.7..locus.tag == locus,]$ACRONYM
            name <- gene
            if(is.na(locus) & length(acronym) == 0) {
                name <- name
            } else if(length(acronym) != 0 & sum(is.na(acronym)) == 0) {
                name <- acronym
            } else if (!is.na(locus)) {
                name <- locus 
            } else {
                name <- name 
            }
            
            par(xpd = F, cex = 0.8)
            barplot(sub_lfc,
                    beside = TRUE,
                    legend = rownames(sub_lfc),
                    xlab = "Contrast",
                    ylab = "LFC",
                    main = paste("LFC of", name, "Isoforms \n", locus),
                    args.legend = list(x = "right",
                                       inset = c(-0.8, 0))
            )
            
            abline(h = 0, col = "black", lwd = 3)
            abline(h = -1, col = "red", lwd = 1)
            abline(h = -0.5, col = "blue", lwd = 1)
            abline(h = 1, col = "red", lwd = 1)
            abline(h = 0.5, col = "blue", lwd = 1)
            abline(v = (length(rownames(sub_lfc)) + 1.5), col = "black", lwd = 1, lty = 2)
            axis(2, 0.5, "0.5", tick = F)
            axis(2, 1, "1", tick = F)
        }
    }
}

convertIsoformNames <- function(conversion_table, isoforms, gene) {
    c_sub <- conversion_table[conversion_table$GENEID == gene,]
    rnames <- rep(NA, length(isoforms))
    count <- 1
    
    if(nrow(c_sub) == 0) {
        return(NA)
    } else {
        for(i in isoforms) {
            for(j in 1:nrow(c_sub)) {
                comp <- c_sub$BAMBUTX[j]
                if(i == comp) {
                    if(is.na(c_sub$acronym[1])) {
                        rnames[count] <- c_sub$ISOFORM[j]
                        count <- count + 1
                    } else {
                        rnames[count] <- c_sub$acronym[j]
                        count <- count + 1
                    }
                }
            }
        }
        return(rnames)
    }
}

calculateReplicateDiffs <- function(cts, taxa) {
    order <- colnames(cts)
    diff_order <- c(
        order[1],
        paste(order[2], "-", order[3]),
        paste(order[2], "-", order[4]),
        paste(order[3], "-", order[4])
    )
    
    diffs <- foreach(
        tx = taxa,
        .combine = "rbind"
    ) %do% {
        txrow <- cts[cts$TXNAME == tx,]
        df <- data.frame(matrix(NA, ncol = ncol(cts), nrow = 1))
        df[1,1] = tx
        df[1,2] = abs(txrow[2] - txrow[3])
        df[1,3] = abs(txrow[2] - txrow[4])
        df[1,4] = abs(txrow[3] - txrow[4])
        colnames(df) <- diff_order 
        df
    }
    return(diffs)
}

getReplicateDiffStats <- function(diffs, plots_only = F) {
    title <- paste("Average Differences in Replicates by Taxa\n ", colnames(diffs)[2], "|", colnames(diffs)[3], "|", colnames(diffs)[4])
    xl <- "Taxa index"
    yl <- "Average difference in replicates"
    if (is.na(plots_only)) {
        for(i in 1:nrow(diffs)) {
            avg <- mean(as.numeric(diffs[i,2:4]))
            diffs$avg[i] <- avg
        }
        
        dens <- density(as.numeric(diffs$avg))
        return(diffs)
        
    } else if(plots_only == T) {
        dens <- density(as.numeric(diffs$avg))
        hist(diffs$avg, main = title)
        barplot(diffs$avg, xlab = xl, ylab = yl, main = title)
        plot(diffs$avg, pch = 16, xlab = xl, ylab = yl, main = title, cex = 0.5, log = "y")
        abline(h = 10000, lty = 2)
        abline(h = -10000, lty = 2)
        abline(h = 1000, col = "blue")
        abline(h = -1000, col = "blue")
        abline(h = 2000, col = "red")
        abline(h = -2000, col = "red")
        plot(dens, main = paste("Density of", title))
        return(diffs)
    } else {
        for(i in 1:nrow(diffs)) {
            avg <- mean(as.numeric(diffs[i,2:4]))
            diffs$avg[i] <- avg
        }
        
        dens <- density(as.numeric(diffs$avg))
        hist(diffs$avg)
        barplot(diffs$avg, xlab = xl, ylab = yl, main = title)
        plot(diffs$avg, pch = 16, xlab = xl, ylab = yl, main = title, cex = 0.5, log = "y")
        abline(h = 10000, lty = 2)
        abline(h = -10000, lty = 2)
        abline(h = 1000, col = "blue")
        abline(h = -1000, col = "blue")
        abline(h = 2000, col = "red")
        abline(h = -2000, col = "red")
        plot(dens)
        return(diffs)
    }
}       

plotReplicateDiagnostics <- function(diff_table, c_table, threshold) {
    extremes <- diff_table[diff_table$avg > threshold,] 
    
    extreme_genes <- c()
    count <- 1
    for(tx in extremes$TXNAME) {
        extreme_genes[count] <- c_table[c_table$BAMBUTX == tx,]$GENEID
        count <- count + 1
    }
    extreme_genes <- unique(extreme_genes)
    
    plot(NA, xlim = c(-7000, 250000), ylim = c(0, 0.0004),
         main = paste("Density of Differences for each Taxa in\n", colnames(diff_table)[2], "|", colnames(diff_table)[3], "|", colnames(diff_table)[4]),
         xlab = "Difference",
         ylab = "Density")
    for(i in 1:nrow(extremes)) {
        lines(density(as.numeric(extremes[i,2:4])))
    }
    lines(density(as.numeric(extremes$avg)), col = "red", lwd = 3)
    abline(v = 1000, col = "darkgreen", lwd = 2)
    abline(v = 2000, col = "blue", lwd = 2)
    abline(v = 5000, col = "darkorange", lwd = 2)
    abline(v = 10000, col = "red", lwd = 2)
    abline(v = 100000, col = "black", lty = 3, lwd = 2)
}

plotRDT <- function(diff_table) {
    rtd <- data.frame(Difference = seq(1000, 10000, by = 1000), Quantity = rep(NA, 10))
    count <- 1
    for(i in rtd$Difference) {
        rtd$Quantity[count] <- nrow(diff_table[diff_table$avg > i,])
        count <- count + 1
    }
    plot(rtd$Difference, rtd$Quantity, pch = 16)
    lines(rtd$Difference, rtd$Quantity, lty = 3)
}

getSharedDE <- function(results, group1, group2) {
    temp <- cbind(rownames(results), results)
    colnames(temp)[1] <- "ID"
    shared <- copyDFscructure(temp, all = F)
    count <- 1
    group1 <- group1 + 1
    group2 <- group2 + 1
    for(i in 1:nrow(results)) {
        if(temp[i, group1] != 0 && temp[i, group2] != 0) {
            shared[count,] <- temp[i,]
            count <- count + 1
        } else {
            next
        }
    }
    rownames(shared) <- shared$ID
    shared <- shared[,2:ncol(shared)]
    return(shared)
}

findGeneFamilies <- function(acro_table, gene_families) {
    big_names <- copyDFscructure(acro_table, all = F) 
    gene_family_counts <- NULL
    gene_count <- NULL
    count <- 1
    
    for(family in gene_families) {
        current_fam <- findGenes(acro_table, family)
        if (!is.na(current_fam[1,1])) {
            gene_count[count] <- nrow(current_fam)
            big_names <- rbind(current_fam, big_names)
            count <- count + 1
        } else {
            gene_count[count] <- 0
            big_names <- rbind(current_fam, big_names)
            count <- count + 1
        }
    }
    
    
    big_names_summary <- data.frame(gene_families,
                                    gene_count)
    
    colnames(big_names_summary) <- c("Gene Family",
                                     "Count"
    )
    
    data <- NULL
    data$big_names <- big_names
    data$big_names_summary <- big_names_summary
    return(data)
}

findGeneFamilies2 <- function(acro_table, gene_families, conversion_table) {
    big_names <- copyDFscructure(acro_table, all = F) 
    gene_family_counts <- NULL
    gene_count <- NULL
    count <- 1
    
    for(family in gene_families) {
        current_fam <- findGenes(acro_table, family)
        if (!is.na(current_fam[1,1])) {
            gene_count[count] <- nrow(current_fam)
            big_names <- rbind(current_fam, big_names)
            count <- count + 1
        } else {
            gene_count[count] <- 0
            big_names <- rbind(current_fam, big_names)
            count <- count + 1
        }
    }
    
    
    
    big_names_summary <- data.frame(gene_families,
                                    gene_count)
    
    colnames(big_names_summary) <- c("Gene Family",
                                     "Count"
    )
    
    data <- NULL
    data$big_names <- na.omit(big_names)
    data$big_names_summary <- big_names_summary
    
    return(data)
}


getDEdiff <- function(gene_set1, gene_set2, gene_fam) {
    all_diffs <- list() 
    all_count <- 1
    
    
    for(gene in gene_fam) {
        diff1 <- NULL
        diff2 <- NULL
        
        set1 <- gene_set1$big_names[grepl(gene, gene_set1$big_names$ACRONYM, ignore.case = T),]
        set2 <- gene_set2$big_names[grepl(gene, gene_set2$big_names$ACRONYM, ignore.case = T),]
        
        count <- 1
        for(i in set1$ACRONYM) {
            if(sum(i == set2) > 0) {
                next
            } else {
                diff1[count] <- i
                count <- count + 1
            }
        }
        
        count <- 1
        for(i in set2$ACRONYM) {
            if(sum(i == set1) > 0) {
                next
            } else {
                diff2[count] <- i
                count <- count + 1
            }
        }
        diffs <- list(diff1, diff2)
        names(diffs) <- c(gene_set1$contrast, gene_set2$contrast)
        all_diffs[[gene]] <- diffs
        all_count <- all_count + 1
    }
    
    names(all_diffs) <- gene_fam
    
    return(all_diffs)
}

getDEgenes <- function(de_table, acrotable, gene_family) {
    all_degenes <- NULL
    for(fam in gene_family) {
        geneset <- acrotable[grepl(fam, acrotable$ACRONYM),][,1]
        geneacro <- acrotable[grepl(fam, acrotable$ACRONYM),][,2]
        
        count <- 1
        for(i in geneset) {
            geneset[count] <- paste("gene_biotype mRNA;", i, sep = " ")
            count <- count + 1
        }
        if(length(geneset) == 0) {
            message("No genes found")
            return(NA)
        }
        degenes <- copyDFscructure(de_table, all = F) 
        count <- 1
        for(i in geneset) {
            rows <- de_table[rownames(de_table) == i,]
            if(nrow(rows) != 0) {
                degenes[count,] <- rows
                count <- count + 1
            } else {
                next
            }
        }
        rownames(degenes) <- geneset 
        degenes$acronym <- geneacro
        all_degenes <- rbind(all_degenes, degenes)
    }
    
    return(all_degenes)
}

getDTUgenes <- function(dtu_table, acrotable, ctable, gene_family) {
    all_degenes <- NULL
    for(fam in gene_family) {
        geneset <- acrotable[grepl(fam, acrotable$ACRONYM),][,1]
        
        count <- 1
        for(i in geneset) {
            geneset[count] <- paste("gene_biotype mRNA;", i, sep = " ")
            count <- count + 1
        }
        
        if(length(geneset) == 0) {
            message("No genes found")
            return(NA)
        }
        txset <- NULL
        txacro <- NULL
        for(gene in geneset) {
            subset <- ctable[ctable$GENEID == gene,][,2:4]
            if(nrow(subset) == 0) {
                message("Gene is not in conversion table")
                return(NA)
            }
            txset <- append(txset, subset$BAMBUTX)
            txacro <- append(txacro, subset$acronym)
        }
        
        txacronym <- NULL
        rnames <- NULL
        degenes <- copyDFscructure(dtu_table, all = F) 
        count <- 1
        acount <- 1
        for(i in txset) {
            rows <- dtu_table[rownames(dtu_table) == i,]
            if(nrow(rows) != 0) {
                degenes[count,] <- rows
                txacronym[count] <- txacro[acount]
                rnames[count] <- i
                count <- count + 1
                acount <- acount + 1
            } else {
                acount <- acount + 1
                next
            }
        }
        rownames(degenes) <- rnames 
        degenes$acronym <- txacronym
        all_degenes <- rbind(all_degenes, degenes)
    }
    
    return(all_degenes)
}

b2n <- function(ids) {
    new_ids <- NULL
    count <- 1
    for(i in ids) {
        split <- strsplit(i, ";")[[1]]
        if(length(split) == 1) {
            new_ids[count] <- trimws(split[1])
            count <- count + 1
        } else {
            new_ids[count] <- trimws(split[2])
            count <- count + 1
        }
    }
    return(new_ids)
}

n2b <- function(ids) {
    new_ids <- NULL
    count <- 1
    for(i in ids) {
        new_ids[count] <- paste()
    }
}

id2acro <- function(ids, acrotable) {
    acros <- NULL
    count <- 1
    for(id in ids) {
        acros[count] <- acrotable[acrotable[,1] == id,2]
        count <- count + 1
    }
    return(acros)
}

plotLFC <- function(deLfc, degenes, title = NULL, cols = NULL, minlfc = 0, acrotable = NULL, hide_legend = F) {
    lfcGene <- copyDFscructure(deLfc)
    rnames <- NULL
    count <- 1
    for(id in rownames(degenes)) {
        curow <- deLfc[rownames(deLfc) == id,]
        if(nrow(curow) != 0) {
            lfcGene[count,] <- curow
            rnames[count] <- rownames(curow)
            count <- count + 1
        } else {
            next
        }
    }
    rownames(lfcGene) <- rnames
    # print(lfcGene)
    # max <- 0 + minlfc
    # min <- 0 - minlfc
    # lfcGene <- lfcGene[lfcGene > min & lfcGene < max]
    # print(lfcGene)
    
    if (!is.null(cols)) {
        bpcols <- cols
    } else {
        bpcols <- c()
        for(i in 1:nrow(lfcGene)) {
            bpcols <- append(bpcols, colors()[i*4])
        }
    }
    par(mai = c(1,1,0.75,2), xpd = T)
    barplot(as.matrix(lfcGene),
            col = bpcols,
            beside = T,
            names.arg = contrasts,
            cex.names = 0.8,
            xlab = "Contrasts",
            ylab = "log2 fold change (LFC)",
            main = title
    )
    
    legendnames <- rownames(lfcGene)
    if(!is.null(acrotable)) {
        legendnames <- id2acro(b2n(rownames(lfcGene)), acrotable)
    }
    textcex = 0.8
    if(hide_legend == F) {
        if(length(legendnames) > 20) {
            textcex <- 0.6
        }
        legend("topright",
               col = bpcols,
               legend = legendnames,
               pch = 15,
               inset = c(-0.31,0),
               cex = textcex
        )
    }
}

getSequenceContent <- function(sequence, type = "DNA") {
    if(type == "DNA") {
        contents <- list(
            "A" = 0,
            "T" = 0,
            "C" = 0,
            "G" = 0,
            total = length(sequence[[1]])
        )
        for(i in 1:length(sequence[[1]])) {
            switch(as.character(sequence[[1]][i]),
                   "A" = contents$A <- contents$A + 1,
                   "T" = contents$T <- contents$T + 1,
                   "C" = contents$C <- contents$C + 1,
                   "G" = contents$G <- contents$G + 1,
            )
        }
        
        return(contents)
    } else if(type == "AA") {
        contents <- list(
            "R" = 0,
            "H" = 0,
            "K" = 0, 
            "D" = 0,
            "E" = 0,
            "S" = 0,
            "T" = 0,
            "N" = 0,
            "Q" = 0,
            "C" = 0,
            "G" = 0,
            "P" = 0,
            "A" = 0,
            "V" = 0,
            "I" = 0,
            "L" = 0,
            "M" = 0,
            "F" = 0,
            "Y" = 0,
            "W" = 0,
            total = length(sequence[[1]])
        )
        for(i in 1:length(sequence[[1]])) {
            switch(as.character(sequence[[1]][i]),
                   "R" = contents$R <- contents$R + 1,
                   "H" = contents$H <- contents$H + 1,
                   "K" = contents$K <- contents$K + 1,
                   "D" = contents$D <- contents$D + 1,
                   "E" = contents$E <- contents$E + 1,
                   "S" = contents$S <- contents$S + 1,
                   "T" = contents$T <- contents$T + 1,
                   "N" = contents$N <- contents$N + 1,
                   "Q" = contents$Q <- contents$Q + 1,
                   "C" = contents$C <- contents$C + 1,
                   "G" = contents$G <- contents$G + 1,
                   "P" = contents$P <- contents$P + 1,
                   "A" = contents$A <- contents$A + 1,
                   "V" = contents$V <- contents$V + 1,
                   "I" = contents$I <- contents$I + 1,
                   "L" = contents$L <- contents$L + 1,
                   "M" = contents$M <- contents$M + 1,
                   "F" = contents$F <- contents$F + 1,
                   "Y" = contents$Y <- contents$Y + 1,
                   "W" = contents$W <- contents$W + 1
            )
        }
        
        return(contents)
    }
}

getSequencePercent <- function(AAcontents, letter) {
    numer <- eval(parse(text = paste(deparse(substitute(AAcontents)), "$", letter, sep = "")))
    return((numer/AAcontents$total) * 100)
}
# These are utility functions
library(RColorBrewer)
library(foreach)

plotSpliceReg <- function(data, set, GeneID, order = NULL) {
    par(mai = c(1.02,0.82,0.82,0.42), xpd = F)
    pre_isoforms <- data[data$GENEID == GeneID,]
    print(pre_isoforms)
    if (!is.null(order)) {
        t_list <- rev((strsplit(trimws(order), " "))[[1]])
        isoforms <- data.frame(matrix(NA, nrow = 1, ncol = ncol(pre_isoforms))) 
        n <- 0
        for(ts in t_list) {
            print(ts)
            n <- n + 1
            isoforms[n,] <- pre_isoforms[rownames(pre_isoforms) == ts,]
        }
        rownames(isoforms) <- isoforms[,1]
        colnames(isoforms) <- colnames(pre_isoforms)
    } else {
        isoforms <- pre_isoforms
    }
    barplot(isoforms$logFC,
            main = paste0("LFC of ", GeneID, "\n", set),
            xlab = "Transcript",
            ylab = "log fold change (LFC)"
            )
    abline(h = 0, col = "red", lwd = 3)
    if (nrow(isoforms) > 3) {
        font.size <- 0.6
    } else {
        font.size <- 0.8
    }
    for (i in 1:nrow(isoforms)) {
        axis(1, i, rownames(isoforms)[i], las = 1, cex.axis = font.size)
    }
}

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

plotIsoform <- function(gene, annotation, exon_marker = F, prop = NULL, acronym_list = NULL, suppress_plot = F) {
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
             labels = transcripts,
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
             labels = transcripts,
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
        print(plist)
        count <- 1
        #fix ordering of props
        for (i in transcripts) {
            for (j in 1:nrow(plist)) {
                if (i == plist$TXNAME[j]) {
                    props[count] <- round(as.double(plist$prop[j]), digit = 4)
                    print(props[count])
                    count <- count + 1
                } else {
                    next
                }
            }
        }
        if (length(props) > 4) {
            if (length(props) >= 6) {
                font.size = 0.6
            } else {
                font.size = 0.8
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
              adj = 1.3,
              padj = -1.2
        )
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
                    strand = featureFrame$strand[1],
                    props = plist
                    )
        )
    } else { 
        transcripts <- unique(featureFrame$transcriptid)
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

plotHeatmapIso <- function(propTable, gene, color) {
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
        colnames(propMat) <- colnames(subset[,3:ncol(subset)])
        
        divx <- seq(1, ncol(propMat))
        divy <- seq(1, nrow(propMat))
        xlimit = c(0, ncol(propMat))
        ylimit = c(0, nrow(propMat))
        par(xpd = T, mai = c(0.5, 2, 0.3, 1))
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
        nums <- round(seq(0, max(propMat), length.out = 10), digit = 2)
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
        # legend(x = xleg,
        #        y = yleg,
        #        inset = c(-0.2,0),
        #        legend = nums,
        #        col = lcol,
        #        pch = 15,
        #        pt.cex = 3
        # )
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
                print(row[i])
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
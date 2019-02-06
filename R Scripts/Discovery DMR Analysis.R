# Discovery DMR Analysis ####
# Diagnosis and Sex
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/2/19

# Packages ####
sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
         "annotatr", "DMRichR"), require, character.only = TRUE)

# Functions ####
CMplotDMR <- function(candidates, prefix, plot.type, bin.max = 450){
        candidates <- cbind("region" = "region", candidates[,c("chr", "start", "pval")], stringsAsFactors = FALSE)
        candidates$chr <- substring(candidates$chr, 4)
        colnames(candidates)[colnames(candidates) == "pval"] <- "pvalue"
        if(plot.type == "m"){
                message("Creating Manhattan plot")
                pdf(paste(prefix, "Manhattan Plot.pdf", sep = " "), height = 5, width = 12)
                par(mar = c(5, 5.5, 1, 3.5), xaxs = "i")
                CMplot(candidates, col = suppressMessages(gg_color_hue(2)), bin.size = 1e7, bin.max = bin.max, cex.axis = 1.2, 
                       plot.type = "m", threshold = 0.05, threshold.lwd = 2, threshold.col = "black", cex = 0.3, 
                       amplify = FALSE, chr.den.col = brewer.pal(9, "YlOrRd"), file.output = FALSE, verbose = FALSE)
                dev.off()
        }
        else{
                if(plot.type == "q"){
                        message("Creating QQ plot")
                        pdf(paste(prefix, "QQ Plot.pdf", sep = " "), height = 5.5, width = 5.5)
                        par(mar = c(4, 4.5, 2.5, 1.5), xaxs = "i", yaxs = "i", mgp = c(2.5, 1, 0))
                        CMplot(manhattan, col = "black", cex.axis = 1.2, plot.type = "q", cex = 0.3, file.output = FALSE, 
                               verbose = FALSE, box = TRUE)
                        dev.off()
                }
                else{
                        print("plot.type must be either m or q.")
                }
        }
}

loadRegions <- function(file, chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE){
        if(grepl("txt", file, fixed = TRUE)){
                regions <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        }
        else{
                regions <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
        }
        if("seqnames" %in% colnames(regions)){
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        }
        regions <- subset(regions, chr %in% chroms)
        regions$chr <- factor(regions$chr, levels = chroms)
        if(sort){
                regions <- regions[order(regions$chr, regions$start),]
        }
        return(regions)
}

plotDendro <- function(ddata, row = !col, col = !row, labels = col) {
        # plot a dendrogram
        yrange <- range(ddata$segments$y)
        yd <- yrange[2] - yrange[1]
        nc <- max(nchar(as.character(ddata$labels$label)))
        if(row){
                tangle <- 0 
        } 
        else { 
                tangle <- 90 
        }
        tshow <- col
        p <- ggplot() +
                geom_segment(data = segment(ddata), aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.45) +
                labs(x = NULL, y = NULL) + theme_dendro()
        if(row) {
                p <- p +
                        scale_x_continuous(expand = c(0.5/length(ddata$labels$x), 0)) +
                        coord_flip()
        } else {
                p <- p +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = "black"))
        }
        return(p)
}

plotLegend <-function(a.gplot){
        # from http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

buildHeatmap <- function(x, phenoData, hm.colors = c("#0000FF", "Black", "#FF0000"), hm.values = c(0, 0.5, 1), 
                         hm.low, hm.high, pheno.breaks, pheno.values) {
        if(is.null(colnames(x))){
                colnames(x) <- sprintf("col%s", 1:ncol(x))
        }
        if(is.null(rownames(x))){
                rownames(x) <- sprintf("row%s", 1:nrow(x))
        }
        
        # Dendrograms
        row.hc <- hclust(dist(x), "ward.D")
        col.hc <- hclust(dist(t(x)), "ward.D")
        row.dendro <- dendro_data(as.dendrogram(row.hc), type = "rectangle")
        col.dendro <- dendro_data(as.dendrogram(col.hc), type = "rectangle")
        col.plot <- plotDendro(col.dendro, col = TRUE, labels = FALSE) +
                theme(plot.margin = unit(c(0, -1.8, 0, -2), "lines"), axis.text.x = element_blank())
        row.plot <- plotDendro(row.dendro, row = TRUE, labels = FALSE) +
                theme(plot.margin = unit(c(0, 2, 0, 0), "lines"))
        col.ord <- match(col.dendro$labels$label, colnames(x))
        row.ord <- match(row.dendro$labels$label, rownames(x))
        xx <- x[row.ord, col.ord]
        dimnames(xx) <- NULL
        xx <- melt(xx)
        
        # Heatmap
        center.plot <- ggplot(xx, aes(X2,X1)) + 
                geom_tile(aes(fill = value, color = value)) +
                scale_fill_gradientn(colors = hm.colors, values = hm.values, limits = c(hm.low, hm.high), na.value = "black") +
                scale_color_gradientn(colors = hm.colors, values = hm.values, limits = c(hm.low, hm.high), na.value = "black") +
                labs(x = NULL, y = NULL) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0), breaks = NULL) +
                theme(plot.margin = unit(rep(0, 4), "lines"))
        
        # phenoData
        sample.ord <- match(col.dendro$labels$label, as.character(phenoData$Sample))
        phenoData$Sample <- factor(as.character(phenoData$Sample), levels = as.character(phenoData$Sample)[sample.ord], 
                                   ordered = TRUE)
        phenoData <- melt(phenoData, id.vars = "Sample")
        phenoData$variable <- factor(phenoData$variable, levels = rev(unique(phenoData$variable)), ordered = TRUE)
        phenoData.plot <- ggplot(phenoData, aes(Sample, variable)) +
                geom_tile(aes(fill = value, color = value)) +
                scale_x_discrete(expand = c(0, 0)) +
                scale_y_discrete(expand = c(0, 0)) +
                scale_color_manual(breaks = pheno.breaks, values = pheno.values) +
                scale_fill_manual(breaks = pheno.breaks, values = pheno.values)
        ret <- list(col = col.plot, row = row.plot, center = center.plot, phenoData = phenoData.plot)
        invisible(ret)
}

printHeatmap <- function(L, widths = c(0.02, 0.8, 0.16, 0.02), heights = c(0.02, 0.15, 0.06, 0.75, 0.02),
                         heatmap.legend.position = c(0.44, -0.45), pheno.legend.position = c(0.925, 0.915)){
        grid.newpage()
        top.layout <- grid.layout(nrow = 5, ncol = 4, widths = unit(widths, "null"), heights = unit(heights, "null"))
        pushViewport(viewport(layout = top.layout))
        
        # Dendrograms
        print(L$col, vp = viewport(layout.pos.col = 2, layout.pos.row = 2))
        print(L$row, vp = viewport(layout.pos.col = 3, layout.pos.row = 4))
        
        # Heatmap
        print(L$center +
                      theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
                            axis.title = element_blank(), legend.position = "none", panel.background = element_blank(), 
                            panel.border = element_blank(), panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), plot.background = element_blank()),
              vp = viewport(layout.pos.col = 2, layout.pos.row = 4))
        
        # PhenoData
        print(L$phenoData +
                      theme_bw(base_size = 24) +
                      theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
                            legend.key = element_blank(), legend.key.size = unit(1, "lines"), 
                            panel.grid.minor = element_blank(), legend.position = "none", 
                            legend.background = element_blank(), legend.text = element_text(size = 12, color = "Black"),
                            plot.margin = unit(c(0, 0, 0, -0.45), "lines"), axis.text = element_blank(),
                            axis.ticks = element_blank(), axis.title = element_blank(), legend.title = element_blank(),
                            plot.title = element_blank()), 
              vp = viewport(layout.pos.col = 2, layout.pos.row = 3))
        
        # Heatmap Legend
        legend <- plotLegend(L$center +
                                     theme(legend.title = element_blank(), 
                                           legend.text = element_text(size = 15, hjust=1, 
                                                                      margin = unit(c(0, 0, 0, 0.2), "lines")),
                                           legend.background = element_blank(), legend.position = heatmap.legend.position))
        pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 2))
        grid.draw(legend)
        upViewport(0)
        
        # PhenoData Legend
        phenoLegend <- plotLegend(L$phenoData +
                                          theme(legend.title = element_blank(), 
                                                legend.text = element_text(size = 15, hjust = 0, 
                                                                           margin = unit(c(0, 0, 0, 0.5), "lines")),
                                                legend.direction = "vertical", legend.position = pheno.legend.position,
                                                legend.background = element_blank()))
        pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 3))
        grid.draw(phenoLegend)
        upViewport(0)
}

ggbiplotPCA <- function(data.pca, groups, pc = 1:2, file, xlim = NULL, ylim = NULL, breaks = c("TD", "ASD"), 
                        values = c("TD" = "#3366CC", "ASD" = "#FF3366"), legend.position = c(0.86, 1.03)){
        # Plots principal components colored by grouping variable and writes file
        pc1 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[1], sep = "")] * 100
        pc2 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[2], sep = "")] * 100
        g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = groups, ellipse = TRUE, circle = FALSE, 
                      var.axes = FALSE, varname.abbrev = FALSE, choices = pc, ellipse.prob = 0.95)
        suppressMessages(g + 
                theme_bw(base_size = 25) +
                theme(legend.direction = 'horizontal', legend.position = legend.position, panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), 
                      legend.spacing.x = unit(0.5, "lines"), plot.margin = unit(c(2, 1, 1, 1), "lines")) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                xlab(paste("PC", pc[1], " (", round(pc1, digits = 0), "% of Variance)", sep = "")) +
                ylab(paste("PC", pc[2], " (", round(pc2, digits = 0), "% of Variance)", sep = "")) +
                scale_color_manual(breaks = breaks , values = values) +
                scale_x_continuous(breaks = pretty_breaks(6)) +
                scale_y_continuous(breaks = pretty_breaks(6)) +
                geom_point(aes(color = groups), size = 3))
        ggsave(file, dpi = 600, width = 8, height = 8, units = "in")
}

methLm <- function(catVars, contVars, sampleData, meth){
        # Analyzes categorical and continuous variables for association with methylation, using linear regression
        stats <- NULL
        methDataCol <- as.numeric(sampleData[,meth])
        for(i in 1:length(catVars)){
                sampleDataCol <- sampleData[,catVars[i]]
                temp <- summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,]
                if(length(temp) == 4){
                        temp <- c(catVars[i], levels(sampleDataCol)[2], temp)
                } 
                else {
                        temp <- cbind(rep(catVars[i], nrow(temp)), 
                                      gsub("sampleDataCol", replacement = "", x = rownames(temp), fixed = TRUE), temp)
                }
                stats <- rbind(stats, temp)
        }
        for(i in 1:length(contVars)){
                sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                temp <- summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,]
                temp <- c(contVars[i], contVars[i], temp)
                stats <- rbind(stats, temp)
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats$Region <- meth
        return(stats)
}

DMRmethLm <- function(DMRs, catVars, contVars, sampleData, file){
        message("Getting DMR methylation by covariate stats")
        start_time <- Sys.time()
        covStats <- lapply(DMRs, function(x){
                methLm(catVars = catVars, contVars = contVars, sampleData = sampleData, meth = x) 
        }) %>% list.rbind
        covStats$Variable <- factor(covStats$Variable, levels = unique(covStats$Variable))
        covStats$Region <- factor(covStats$Region, levels = unique(covStats$Region))
        covStats[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(covStats[,c("Estimate", "StdError", "tvalue", 
                                                                                       "pvalue")], as.numeric)
        covStats$log_pvalue <- -log10(covStats$pvalue)
        covStats$qvalue <- p.adjust(covStats$pvalue, method = "fdr")
        covStats <- covStats[,c("Region", "Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue", "log_pvalue", 
                                "qvalue")]
        message("Complete, writing file")
        write.table(covStats, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
        end_time <- Sys.time()
        message("Time difference of ", round(end_time - start_time, 2), " minutes")
        return(covStats)
}

DMRmethLmSum <- function(covStats, file){
        covSum <- cbind(table(covStats$Variable, covStats$pvalue < 0.05), 
                        table(covStats$Variable, covStats$pvalue < 0.05 & covStats$Estimate > 0)[,"TRUE"],
                        table(covStats$Variable, covStats$pvalue < 0.05 & covStats$Estimate < 0)[,"TRUE"],
                        table(covStats$Variable, covStats$qvalue < 0.05)[,"TRUE"]) %>% as.data.frame
        covSum$Variable <- rownames(covSum)
        rownames(covSum) <- 1:nrow(covSum)
        colnames(covSum) <- c("NoAssoc", "NominalSig", "NominalPos", "NominalNeg", "FDRsig", "Variable")
        covSum$perNominalSig <- covSum$NominalSig * 100 / length(unique(covStats$Region))
        covSum <- covSum[order(-covSum$NominalSig), c("Variable", "NoAssoc", "NominalSig", "perNominalSig", "NominalPos", "NominalNeg", "FDRsig")]
        write.table(covSum, file = file, sep = "\t",
                    quote = FALSE, row.names = FALSE)
        return(covSum)
}

covHeatmap <- function(covStats, variableOrdering = c("unsorted", "manual", "hierarchical"), 
                       regionOrdering = c("unsorted", "variable", "hierarchical"), variables = NULL,
                       sortVariable = "Diagnosis_Alg", file){
        # Sort Variables
        if(variableOrdering == "unsorted"){
                variableOrder <- 1:length(unique(covStats$Variable))
        }
        if(variableOrdering == "manual"){
                variableOrder <- match(variables, unique(covStats$Variable))
        }
        if(variableOrdering == "hierarchical"){
                pvals <- cast(covStats[,c("Region", "Variable", "log_pvalue")], formula = Variable ~ Region, 
                              fun.aggregate = mean, value = "log_pvalue")
                variableOrder <- hclust(dist(pvals[,2:ncol(pvals)], method = "euclidean"), method = "ward.D")$order
        }
        covStats$Variable <- factor(covStats$Variable, levels = unique(covStats$Variable)[variableOrder], 
                                    ordered = TRUE)
        # Sort Regions
        if(regionOrdering == "unsorted"){
                regionOrder <- 1:length(unique(covStats$Region)) 
        }
        if(regionOrdering == "variable"){
                regionOrder <- order(covStats$log_pvalue[covStats$Variable == sortVariable])
        }
        if(regionOrdering == "hierarchical"){
                pvals <- cast(covStats[,c("Region", "Variable", "log_pvalue")], formula = Region ~ Variable, 
                              fun.aggregate = mean, value = "log_pvalue")
                regionOrder <- hclust(dist(pvals[,2:ncol(pvals)], method = "euclidean"), method = "ward.D")$order
        }
        covStats$Region <- factor(covStats$Region, levels = unique(covStats$Region)[regionOrder], 
                                  ordered = TRUE)
        # Plot Heatmap
        gg <- ggplot(data = covStats)
        gg +
                geom_tile(aes(y = Region, x = Variable, fill = log_pvalue)) +
                scale_fill_gradientn("-log(p-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                                     limits = c(0,quantile(x = covStats$log_pvalue, probs = 0.999, names = FALSE))) +
                theme_bw(base_size = 24) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      plot.margin = unit(c(1, 6, 1, 1), "lines"), axis.ticks.x = element_line(size = 1.25), 
                      axis.ticks.y = element_blank(), 
                      axis.text.x = element_text(size = 9, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
                      axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
                      legend.position = c(1.08, 0.84), legend.background = element_blank(), 
                      legend.key.size = unit(1, "lines"), legend.title = element_text(size = 13), 
                      legend.text = element_text(size = 11))
        ggsave(file, dpi = 600, width = 9, height = 6, units = "in")
}

# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_smoothed_methylation_Dx_Discovery50_females.txt",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
meth <- DMRs[,c("chr", "start", "end", "DMRid", colnames(DMRs)[grepl("JLCM", colnames(DMRs), fixed = TRUE)])]
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Candidates
candidates <- loadRegions("DMRs/Discovery/Diagnosis Females 50/CandidateRegions_Dx_Discovery50_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(meth))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis Females DMRs", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis Females DMRs", plot.type = "q")

# Heatmap ####
# Meth Data
meth_heat <- meth[,grepl("JLCM", colnames(meth), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.99, names = FALSE) %>% ceiling # 25

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"), ordered = TRUE)

# Plot Heatmap
pdf(file="Figures/Diagnosis Females DMRs Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData, hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# PCA ####
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
data.pca <- prcomp(meth_heat %>% as.matrix %>% t, center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Females Diagnosis DMRs PCA by Diagnosis.png", xlim = c(-40, 40), ylim = c(-40, 40))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Females Diagnosis DMRs PCA by Study.png", xlim = c(-40, 40), ylim = c(-40, 40),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Covariate Association ####
# Prep Data
meth_cov <- meth_heat %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Sex")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex")]

# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Females Diagnosis DMRs Methylation by Covariate Stats.txt")
covSum <- DMRmethLmSum(covStats, file = "Tables/Females Diagnosis DMRs Methylation by Covariate Summary.txt")
covSum$Variable[covSum$FDRsig > 0]
# [1] "percent_cpg_meth"  "percent_chh_meth"  "GDM"               "C_coverage"        "dedup_reads_M"    
# [6] "CG_coverage"       "percent_duplicate" "Site" 

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Study", "Site", "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", 
             "MSLrlTscore36", "MSLvrTscore36", "ga_w", "bw_g", "percent_trimmed", "percent_aligned",
             "percent_duplicate", "dedup_reads_M", "C_coverage", "CG_coverage", "percent_cpg_meth", "percent_chg_meth",
             "percent_chh_meth", "MomEdu_detail", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", 
             "Mat_Weight_kg_PrePreg", "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy",
             "cotinine_urine_ngml", "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", 
             "AllEQ_tot_All_FA_mcg_Mo_3", "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", 
             "AllEQ_tot_All_FA_mcg_Mo1", "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", 
             "AllEQ_tot_All_FA_mcg_Mo4", "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", 
             "AllEQ_tot_All_FA_mcg_Mo7", "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")

covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Females Diagnosis DMRs Covariate Heatmap Sorted by Diagnosis.png")

covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Females Diagnosis DMRs Covariate Heatmap Clustered.png")

# Gene Annotation ####
# GREAT style
# Also see DM.R, annotatr


# CpG Annotation ####
# See DM.R, annotatr


# GREAT Analysis ####






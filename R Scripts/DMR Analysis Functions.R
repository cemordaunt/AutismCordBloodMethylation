# DMR Analysis Functions ####
# Charles Mordaunt
# 6/7/19

#sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
#         "annotatr", "GenomicRanges", "LOLA", "rtracklayer", "R.utils", "rGREAT", "DMRichR"), require, character.only = TRUE)

# Region Helper Functions -------------------------------------------------

loadRegions <- function(file, chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE, 
                        DMRid = FALSE, DMBid = FALSE){
        if(grepl("txt", file, fixed = TRUE)){
                regions <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        }
        else{
                regions <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
        }
        if(nrow(regions) == 0){
                stop(c(paste("There are no regions in", file, sep = " ")))
        }
        if("seqnames" %in% colnames(regions)){
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        }
        regions <- subset(regions, chr %in% chroms)
        regions$chr <- factor(regions$chr, levels = chroms)
        if(sort){
                regions <- regions[order(regions$chr, regions$start),]
        }
        if(DMRid){
                regions$DMRid <- paste("DMR", 1:nrow(regions), sep = "_")
        }
        if(DMBid){
                regions$DMBid <- paste("DMB", 1:nrow(regions), sep = "_")
        }
        return(regions)
}

writeBED <- function(regions, file){
        if("seqnames" %in% colnames(regions)){
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        }
        if("DMRid" %in% colnames(regions)){
                bed <- data.frame("chr" = regions$chr, "start" = regions$start, "end" = regions$end, "name" = regions$DMRid, 
                                  "score" = 0, "strand" = 0, "thickStart" = 0, "thickEnd" = 0, 
                                  "RGB" = ifelse(regions$percentDifference > 0, "255,0,0", "0,0,255"))
        } else {
                if("DMBid" %in% colnames(regions)){
                        bed <- data.frame("chr" = regions$chr, "start" = regions$start, "end" = regions$end, "name" = regions$DMBid, 
                                          "score" = 0, "strand" = 0, "thickStart" = 0, "thickEnd" = 0, 
                                          "RGB" = ifelse(regions$percentDifference > 0, "255,0,0", "0,0,255"))
                } else {
                        if("percentDifference" %in% colnames(regions)){
                                bed <- data.frame("chr" = regions$chr, "start" = regions$start, "end" = regions$end, 
                                                  "name" = paste("region", 1:nrow(regions), sep = "_"), 
                                                  "score" = 0, "strand" = 0, "thickStart" = 0, "thickEnd" = 0, 
                                                  "RGB" = ifelse(regions$percentDifference > 0, "255,0,0", "0,0,255"))
                        } else {
                                bed <- data.frame("chr" = regions$chr, "start" = regions$start, "end" = regions$end, 
                                                  "name" = paste("region", 1:nrow(regions), sep = "_"), 
                                                  "score" = 0, "strand" = 0, "thickStart" = 0, "thickEnd" = 0)
                        }
                }
        }
        write.table(bed, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

makeGRange <- function(DMRs, direction = c("all", "hyper", "hypo")){
        if(direction == "hyper"){
                DMRs <- subset(DMRs, percentDifference > 0)
        }
        if(direction == "hypo"){
                DMRs <- subset(DMRs, percentDifference < 0)
        }
        if("chr" %in% colnames(DMRs)){
                GR <- GRanges(seqnames = DMRs$chr, ranges = IRanges(start = DMRs$start, end = DMRs$end))
        } else {
                GR <- GRanges(seqnames = DMRs$seqnames, ranges = IRanges(start = DMRs$start, end = DMRs$end))
        }
        return(GR)
}

GRangeExtend <- function(x, extend){
        ranges(x) <- IRanges(start(x) - extend, end(x) + extend)
        x <- trim(x)
        return(x)
}

getRegionStats <- function(DMRs, background, n){
        DMRstats <- sapply(DMRs, function(x){
                c("DMR_Number" = nrow(x), "DMR_Width_KB" = sum(x$width) / 1000, "DMR_CpGs" = sum(x[,colnames(x) %in% c("L", "CpGs")]))
        }) %>% t %>% as.data.frame
        backgroundStats <- sapply(background, function(x){
                c("Background_Number_K" = nrow(x) / 1000, "Background_Width_GB" = sum(x$width) / 1e9, "Background_CpGs_M" = sum(x$n) / 1e6)
        }) %>% t %>% as.data.frame
        regionStats <- cbind(DMRstats, backgroundStats)
        regionStats$Comparison <- row.names(DMRstats)
        row.names(regionStats) <- 1:nrow(regionStats)
        regionStats$n <- n
        regionStats <- regionStats[,c("Comparison", "n", "DMR_Number", "DMR_Width_KB", "DMR_CpGs", "Background_Number_K", 
                                      "Background_Width_GB", "Background_CpGs_M")]
        return(regionStats)
}

# Visualization Functions -------------------------------------------------

CMplotDMR <- function(candidates, prefix, plot.type, bin.max = 450, m.cex = 0.3){
        candidates <- cbind("region" = "region", candidates[,c("chr", "start", "pval")], stringsAsFactors = FALSE)
        candidates$chr <- substring(candidates$chr, 4)
        colnames(candidates)[colnames(candidates) == "pval"] <- "pvalue"
        if(plot.type == "m"){
                message("Creating Manhattan plot")
                pdf(paste(prefix, "Manhattan Plot.pdf", sep = " "), height = 5, width = 12)
                par(mar = c(5, 5.5, 1, 3.5), xaxs = "i")
                CMplot(candidates, col = suppressMessages(gg_color_hue(2)), bin.size = 1e7, bin.max = bin.max, cex.axis = 1.2, 
                       plot.type = "m", threshold = 0.05, threshold.lwd = 2, threshold.col = "black", cex = m.cex, 
                       amplify = FALSE, chr.den.col = brewer.pal(9, "YlOrRd"), file.output = FALSE, verbose = FALSE)
                dev.off()
        }
        else{
                if(plot.type == "q"){
                        message("Creating QQ plot")
                        pdf(paste(prefix, "QQ Plot.pdf", sep = " "), height = 5.5, width = 5.5)
                        par(mar = c(4, 4.5, 2.5, 1.5), xaxs = "i", yaxs = "i", mgp = c(2.5, 1, 0))
                        CMplot(candidates, col = "black", cex.axis = 1.2, plot.type = "q", cex = 0.3, file.output = FALSE, 
                               verbose = FALSE, box = TRUE)
                        dev.off()
                }
                else{
                        print("plot.type must be either m or q.")
                }
        }
}

ggBoxPlot <- function(data, x, y, fill, ylab, legend.name, facet, file, width = 10, height = 7, 
                      legend.position = c(0.85, 0.38), legend.direction = "vertical", 
                      axis.ticks.x = element_blank(), outlier.size = 0.8,
                      axis.text.x = element_blank(), nrow = 2, ncol = NULL, ylim = c(0,100)){
        gg <- ggplot(data = data)
        gg <- gg +
                geom_boxplot(aes(x = x, y = y, fill = fill), size = 0.8, outlier.size = outlier.size) +
                theme_bw(base_size = 24) +
                theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
                      axis.ticks.x = axis.ticks.x, legend.key = element_blank(), panel.grid.minor = element_blank(),
                      legend.position = legend.position, legend.background = element_blank(), axis.text.x = axis.text.x, 
                      legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 19), 
                      axis.ticks.y = element_line(size = 1.25), legend.title = element_text(size = 22),
                      strip.background = element_blank(), legend.direction = legend.direction, panel.spacing.y = unit(0, "lines"), 
                      plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title.x = element_blank(), 
                      axis.text.y = element_text(size = 16, color = "black")) +
                ylab(ylab) +
                scale_fill_manual(name = legend.name, values = c("#3366CC", "#FF3366")) +
                scale_color_manual(name = legend.name, values = c("#3366CC", "#FF3366")) +
                coord_cartesian(ylim = ylim) +
                facet_wrap(facets = facet, nrow = nrow, ncol = ncol)
        ggsave(filename = file, plot = gg, dpi = 600, width = width, height = height, units = "in")
}

ggScatterPlot <- function(x, y, groupVar, fileName, xlab, ylab, xlim = NULL, ylim = NULL, legendPos = c(0.87,1.03)){
        # Plots 2 continuous variables colored by a grouping variable and writes the file
        g <- ggplot()
        g + 
                geom_smooth(aes(x=x, y=y), method="lm") +
                geom_point(aes(x=x, y=y, color=groupVar), size=3) +
                theme_bw(base_size = 25) +
                theme(legend.direction = 'horizontal', legend.position = legendPos, panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), 
                      plot.margin = unit(c(2,1,1,1), "lines")) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                scale_x_continuous(breaks=pretty_breaks(n=5)) +
                scale_y_continuous(breaks=pretty_breaks(n=5)) +
                xlab(xlab) +
                ylab(ylab) +
                scale_color_manual(breaks = c(levels(groupVar)[1], levels(groupVar)[2]), values = c("#3366CC", "#FF3366"))
        ggsave(fileName, dpi = 600, width = 8, height = 7, units = "in")
}

.plotDendro <- function(ddata, row = !col, col = !row, labels = col) {
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

.plotLegend <-function(a.gplot){
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
        col.plot <- .plotDendro(col.dendro, col = TRUE, labels = FALSE) +
                theme(plot.margin = unit(c(0, -1.8, 0, -2), "lines"), axis.text.x = element_blank())
        row.plot <- .plotDendro(row.dendro, row = TRUE, labels = FALSE) +
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
        legend <- .plotLegend(L$center +
                                     theme(legend.title = element_blank(), 
                                           legend.text = element_text(size = 15, hjust=1, 
                                                                      margin = unit(c(0, 0, 0, 0.2), "lines")),
                                           legend.background = element_blank(), legend.position = heatmap.legend.position))
        pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 2))
        grid.draw(legend)
        upViewport(0)
        
        # PhenoData Legend
        phenoLegend <- .plotLegend(L$phenoData +
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

# Covariate Association Functions -----------------------------------------

.methLm <- function(catVars, contVars, sampleData, meth, adj = NULL){
        # Analyzes categorical and continuous variables for association with methylation, 
        # using linear regression with an optional adjustment variable
        stats <- NULL
        methDataCol <- as.numeric(sampleData[,meth])
        if(!is.null(adj)){
                adjDataCol <- sampleData[,adj]
        }
        if(!is.null(catVars)){
                for(i in 1:length(catVars)){
                        sampleDataCol <- sampleData[,catVars[i]]
                        if(!is.null(adj)){
                                temp <- tryCatch(summary(lm(methDataCol ~ sampleDataCol + adjDataCol))$coefficients[-1,],
                                                 error = function(x){
                                                         message("Error with ", meth, " and ", catVars[i])
                                                         rep(NA, 4) # Returns NA if error in lm
                                                 })
                                if(!is.na(temp[1])){temp <- temp[-nrow(temp),]} # Remove adj var if temp is not NA
                        }
                        else {
                                temp <- tryCatch(summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,],
                                                 error = function(x){
                                                         message("Error with ", meth, " and ", catVars[i])
                                                         rep(NA, 4) # Returns NA if error in lm
                                                 })
                        }
                        if(length(temp) == 4){
                                temp <- c(catVars[i], levels(sampleDataCol)[2], temp)
                        } 
                        else {
                                temp <- cbind(rep(catVars[i], nrow(temp)), 
                                              gsub("sampleDataCol", replacement = "", x = rownames(temp), fixed = TRUE), 
                                              temp)
                        }
                        stats <- rbind(stats, temp)
                }
        }
        if(!is.null(contVars)){
                for(i in 1:length(contVars)){
                        sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                        sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                        if(!is.null(adj)){
                                temp <- summary(lm(methDataCol ~ sampleDataCol + adjDataCol))$coefficients[-1,]
                                temp <- temp[-nrow(temp),] # Remove adj var
                        }
                        else {
                                temp <- summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,]
                        }
                        temp <- c(contVars[i], contVars[i], temp)
                        stats <- rbind(stats, temp)
                }
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats$Region <- meth
        return(stats)
}

DMRmethLm <- function(DMRs, catVars, contVars, sampleData, file, adj = NULL){
        if(!is.null(adj)){
                message("[DMRmethLm] Getting DMR methylation by covariate stats using adjustment variable")
        }
        else {
                message("[DMRmethLm] Getting DMR methylation by covariate stats")
        }
        covStats <- lapply(DMRs, function(x){
                .methLm(catVars = catVars, contVars = contVars, sampleData = sampleData, meth = x, adj = adj) 
        }) %>% list.rbind
        covStats$Variable <- as.character(covStats$Variable)
        covStats$Variable[covStats$Variable %in% c("Site", "MomEdu_detail")] <- paste(covStats$Variable[covStats$Variable %in% c("Site", "MomEdu_detail")],
                                                                                      covStats$Term[covStats$Variable %in% c("Site", "MomEdu_detail")], sep = "_")
        covStats$Variable <- gsub(" ", "_", covStats$Variable)
        covStats$Variable <- factor(covStats$Variable, levels = unique(covStats$Variable))
        covStats$Region <- factor(covStats$Region, levels = unique(covStats$Region))
        covStats[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(covStats[,c("Estimate", "StdError", "tvalue", 
                                                                                       "pvalue")], as.numeric)
        covStats$log_pvalue <- -log10(covStats$pvalue)
        covStats$qvalue <- p.adjust(covStats$pvalue, method = "fdr")
        covStats <- covStats[,c("Region", "Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue", "log_pvalue", 
                                "qvalue")]
        message("[DMRmethLm] Complete, writing file")
        write.table(covStats, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
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
                       sortVariable = "Diagnosis_Alg", file, probs = 0.999, axis.text.y = element_blank(), 
                       legend.position = c(1.07, 0.81), width = 9.5, height = 6, axis.text.x.size = 10){
        # Replace NA/Inf/NaN Values with 0
        covStats$log_pvalue[is.na(covStats$log_pvalue) | is.infinite(covStats$log_pvalue) | 
                                    is.nan(covStats$log_pvalue)] <- 0
        # Sort Variables
        if(variableOrdering == "unsorted"){
                variableOrder <- 1:length(unique(covStats$Variable))
        }
        if(variableOrdering == "manual"){
                variableOrder <- match(variables, unique(covStats$Variable))
        }
        if(variableOrdering == "hierarchical"){
                pvals <- reshape::cast(covStats[,c("Region", "Variable", "log_pvalue")], formula = Variable ~ Region, 
                              fun.aggregate = mean, value = "log_pvalue", add.missing = TRUE, fill = 0) # Cast sorts variables
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
                pvals <- reshape::cast(covStats[,c("Region", "Variable", "log_pvalue")], formula = Region ~ Variable, 
                              fun.aggregate = mean, value = "log_pvalue", add.missing = TRUE, fill = 0)
                regionOrder <- hclust(dist(pvals[,2:ncol(pvals)], method = "euclidean"), method = "ward.D")$order
        }
        covStats$Region <- factor(covStats$Region, levels = unique(covStats$Region)[regionOrder], 
                                  ordered = TRUE)
        
        # Label Replacements
        replacements <- c("Site_Drexel" = "Site (Drexel)", "Site_Johns_Hopkins_University" = "Site (Johns Hopkins)", 
                          "Site_Kaiser_Permanente" = "Site (Kaiser)", "Diagnosis_Alg" = "Diagnosis", "ADOScs" = "ADOS", 
                          "MSLelcStandard36" = "Mullen Composite", "MSLelTscore36" = "Mullen Expressive Language", 
                          "MSLfmTscore36" = "Mullen Fine Motor", "MSLrlTscore36" = "Mullen Receptive Language", 
                          "MSLvrTscore36" = "Mullen Visual Reception", "percent_trimmed" = "Bases Trimmed", 
                          "percent_aligned" = "Aligned Reads", "percent_duplicate" = "PCR Duplicates", 
                          "dedup_reads_M" = "Unique Reads", "C_coverage" = "C Coverage", "CG_coverage" = "CpG Coverage", 
                          "percent_chg_meth" = "Global CHG Methylation", "percent_chh_meth" = "Global CHH Methylation", "ga_w" = "Gestational Age",
                          "bw_g" = "Birthweight", "MomAgeYr" = "Maternal Age", "MomEdu_detail_1" = "Maternal Edu (8th Grade)", 
                          "MomEdu_detail_2" = "Maternal Edu (Some High School)", 
                          "MomEdu_detail_3" = "Maternal Edu (High School Diploma)", "MomEdu_detail_4" = "Maternal Edu (Some College)",
                          "MomEdu_detail_5" = "Maternal Edu (Associate's)",  "MomEdu_detail_7" = "Maternal Edu (Master's)", 
                          "MomEdu_detail_8" = "Maternal Edu (Doctorate)", "Mat_Height_cm" = "Maternal Height", 
                          "Mat_Weight_kg_PrePreg" = "Maternal Weight", "Mat_BMI_PrePreg" = "Maternal BMI", 
                          "DM1or2" = "Maternal Diabetes", "GDM" = "Maternal Gestational Diabetes", "PE" = "Maternal Preeclampsia", 
                          "parity" = "Parity", "dad_age" = "Paternal Age", "home_ownership" = "Own Home", 
                          "marital_status" = "Married", "SmokeYN_Pregnancy" = "Maternal Smoking", 
                          "cotinine_urine_ngml" = "Urine Cotinine", "final_creatinine_mgdl" = "Urine Creatinine", 
                          "percent_cpg_meth_bsseq" = "Global CpG Methylation", "Bcell" = "B Cells", "CD4T" = "CD4 T Cells", 
                          "CD8T" = "CD8 T Cells", "Gran" = "Granulocytes", "Mono" = "Monocytes", "NK" = "NK Cells", "nRBC" = "nRBCs")
        
        # Plot Heatmap
        gg <- ggplot(data = covStats)
        gg +
                geom_tile(aes(y = Region, x = Variable, fill = log_pvalue)) +
                scale_fill_gradientn("-log(p-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                                     limits = c(0,quantile(x = covStats$log_pvalue, probs = probs, names = FALSE, na.rm = TRUE))) +
                scale_x_discrete(labels = function(x){str_replace_all(string = x, pattern = replacements)}) +
                theme_bw(base_size = 24) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      plot.margin = unit(c(1, 5.5, 1, 1), "lines"), axis.ticks.x = element_line(size = 1), 
                      axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
                      axis.text.x = element_text(size = axis.text.x.size, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
                      axis.text.y = axis.text.y, axis.title = element_blank(), legend.key = element_blank(),  
                      legend.position = legend.position, legend.background = element_blank(), 
                      legend.key.size = unit(1, "lines"), legend.title = element_text(size = 12), 
                      legend.text = element_text(size = 11))
        ggsave(file, dpi = 600, width = width, height = height, units = "in")
}

# Annotation Functions ----------------------------------------------------

getDMRanno <- function(DMRstats, regDomains, file = NULL){
        message("[getDMRanno] Adding genes to DMRs by regulatory domains")
        GR_regDomains <- GRanges(seqnames = regDomains$gene_chr, ranges = IRanges(start = regDomains$distal_start, end = regDomains$distal_end))
        GR_DMRstats <- GRanges(seqnames = DMRstats$chr, ranges = IRanges(start = DMRstats$start, end = DMRstats$end))
        overlaps <- as.data.frame(findOverlaps(GR_DMRstats, GR_regDomains))
        rm(GR_regDomains, GR_DMRstats)
        if("DMBid" %in% colnames(DMRstats)){
                colnames(DMRstats)[colnames(DMRstats) == "DMBid"] <- "DMRid"
        }
        DMRstats_genes <- cbind("DMRid" = DMRstats$DMRid[overlaps$queryHits], regDomains[overlaps$subjectHits,], row.names=NULL)
        rm(overlaps)
        DMRstats_genes <- merge(DMRstats, DMRstats_genes, by = "DMRid", all = TRUE, sort = FALSE)
        message("[getDMRanno] Getting DMR positions relative to genes")
        message("[getDMRanno] Gene positions added:\t")
        DMRstats_annotated <- NULL
        for(i in 1:nrow(DMRstats_genes)){
                if(i %% 500 == 0){message(i, "\t", appendLF = FALSE)}
                temp <- DMRstats_genes[i,]
                if(!temp$gene_strand %in% c("+", "-")){temp$distanceToTSS <- NA; temp$positionInGene <- NA} 
                else {
                        if(temp$gene_strand == "+"){ # + strand
                                if(temp$start < temp$gene_start & temp$end < temp$gene_start){temp$distanceToTSS <- temp$end - temp$gene_start} #upstream of TSS (-)
                                if(temp$start <= temp$gene_start & temp$end >= temp$gene_start){temp$distanceToTSS <- 0} #overlapping TSS
                                if(temp$start > temp$gene_start & temp$end > temp$gene_start){temp$distanceToTSS <- temp$start - temp$gene_start} #downstream of TSS (+)
                                if(temp$end < temp$gene_start){temp$positionInGene <- "upstream"}
                                if(temp$end > temp$gene_start & temp$start < temp$gene_end){temp$positionInGene <- "gene_body"}
                                if(temp$start > temp$gene_end){temp$positionInGene <- "downstream"}
                                if(temp$distanceToTSS == 0){temp$positionInGene <- "TSS"}
                        }
                        else { # - strand
                                if(temp$start > temp$gene_end & temp$end > temp$gene_end){temp$distanceToTSS <- temp$gene_end - temp$start} #upstream of TSS (-)
                                if(temp$start <= temp$gene_end & temp$end >= temp$gene_end){temp$distanceToTSS <- 0} #overlapping TSS
                                if(temp$start < temp$gene_end & temp$end < temp$gene_end){temp$distanceToTSS <- temp$gene_end - temp$end} #downstream of TSS (+)
                                if(temp$start > temp$gene_end){temp$positionInGene <- "upstream"}
                                if(temp$start < temp$gene_end & temp$end > temp$gene_start){temp$positionInGene <- "gene_body"}
                                if(temp$end < temp$gene_start){temp$positionInGene <- "downstream"}
                                if(temp$distanceToTSS == 0){temp$positionInGene <- "TSS"}
                        }
                }
                DMRstats_annotated <- rbind(DMRstats_annotated, temp)
        }
        DMRstats_annotated$DMRid <- factor(DMRstats_annotated$DMRid, levels = unique(DMRstats_annotated$DMRid))
        DMRstats_annotated <- aggregate(formula = cbind(gene_name, gene_entrezID, gene_strand, distanceToTSS, positionInGene) ~ DMRid, 
                                        data = DMRstats_annotated, FUN = function(x) paste(x, collapse = ", "), simplify = TRUE)
        DMRstats_annotated <- merge(x = DMRstats, y = DMRstats_annotated, by = "DMRid", all.x = TRUE, all.y = FALSE, sort = FALSE)
        message("\n[getDMRanno] Getting CpG annotations")
        annotations <- build_annotations(genome = "hg38", annotations = "hg38_cpgs")
        annotations <- GenomeInfoDb::keepStandardChromosomes(annotations, pruning.mode = "coarse")
        DMRs_CpGs <- annotate_regions(regions = GRanges(seqnames = DMRstats$chr, 
                                                        ranges = IRanges(start = DMRstats$start, end = DMRstats$end), 
                                                        DMRid = DMRstats$DMRid),
                                      annotations = annotations, ignore.strand = TRUE, quiet = TRUE) %>% as.data.frame
        colnames(DMRs_CpGs)[colnames(DMRs_CpGs) == "annot.type"] <- "CpG_Anno"
        DMRs_CpGs$CpG_Anno[DMRs_CpGs$CpG_Anno == "hg38_cpg_islands"] <- "CpG_Island"
        DMRs_CpGs$CpG_Anno[DMRs_CpGs$CpG_Anno == "hg38_cpg_shores"] <- "CpG_Shore"
        DMRs_CpGs$CpG_Anno[DMRs_CpGs$CpG_Anno == "hg38_cpg_shelves"] <- "CpG_Shelf"
        DMRs_CpGs$CpG_Anno[DMRs_CpGs$CpG_Anno == "hg38_cpg_inter"] <- "CpG_Open_Sea"
        DMRs_CpGs <- aggregate(formula = CpG_Anno ~ DMRid, data = DMRs_CpGs, FUN = function(x) paste(unique(x), collapse = ", "), simplify = TRUE)
        DMRstats_annotated <- merge(x = DMRstats_annotated, y = DMRs_CpGs, by = "DMRid", all.x = TRUE, all.y = FALSE, sort = FALSE)
        message("[getDMRanno] Getting gene regulatory annotations")
        annotations <- build_annotations(genome = "hg38", annotations = c("hg38_basicgenes", "hg38_genes_intergenic", 
                                                                          "hg38_genes_intronexonboundaries", 
                                                                          "hg38_enhancers_fantom"))
        annotations <- GenomeInfoDb::keepStandardChromosomes(annotations, pruning.mode = "coarse")
        DMRs_GeneReg <- annotate_regions(regions = GRanges(seqnames = DMRstats$chr, 
                                                           ranges = IRanges(start = DMRstats$start, end = DMRstats$end), 
                                                           DMRid = DMRstats$DMRid),
                                         annotations = annotations, ignore.strand = TRUE, quiet = TRUE) %>% as.data.frame
        colnames(DMRs_GeneReg)[colnames(DMRs_GeneReg) == "annot.type"] <- "GeneReg_Anno"
        DMRs_GeneReg$GeneReg_Anno <- gsub(pattern = "hg38_", replacement = "", x = DMRs_GeneReg$GeneReg_Anno, fixed = TRUE)
        DMRs_GeneReg <- aggregate(formula = GeneReg_Anno ~ DMRid, data = DMRs_GeneReg, FUN = function(x) paste(unique(x), collapse = ", "), simplify = TRUE)
        DMRstats_annotated <- merge(x = DMRstats_annotated, y = DMRs_GeneReg, by = "DMRid", all.x = TRUE, all.y = FALSE, sort = FALSE)
        DMRstats_annotated <- DMRstats_annotated[order(DMRstats_annotated$chr, DMRstats_annotated$start),]
        if(sum(grepl(pattern = "DMB", x = DMRstats_annotated$DMRid, fixed = TRUE)) > 0){
                colnames(DMRstats_annotated)[colnames(DMRstats_annotated) == "DMRid"] <- "DMBid"
        }
        if(!is.null(file)){
                write.table(DMRstats_annotated, file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
        message("[getDMRanno] Complete!")
        return(DMRstats_annotated)
}

getDMRgeneList <- function(DMRstats, regDomains, direction = c("all", "hyper", "hypo"), 
                           type = c("gene_name", "gene_entrezID")){
        if(is.null(direction)){direction <- "all"}
        if(direction == "hyper"){DMRstats <- subset(DMRstats, percentDifference > 0)}
        if(direction == "hypo"){DMRstats <- subset(DMRstats, percentDifference < 0)}
        GR_regDomains <- GRanges(seqnames = regDomains$gene_chr, ranges = IRanges(start = regDomains$distal_start, end = regDomains$distal_end))
        GR_DMRstats <- GRanges(seqnames = DMRstats$chr, ranges = IRanges(start = DMRstats$start, end = DMRstats$end))
        overlaps <- as.data.frame(findOverlaps(GR_DMRstats, GR_regDomains))
        geneList <- regDomains[overlaps$subjectHits, type] %>% unique %>% sort
        message(nrow(DMRstats), " input regions correspond with ", length(geneList), " ", type, "s.")
        return(geneList)
}

# Enrichment Functions ----------------------------------------------------

getDAVID <- function(genes, background, file, categories, benjamini = 0.05){
        message("[getDAVID] Uploading gene and background lists")
        david <- DAVIDWebService$new(url = "https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/",
                                     email = "cemordaunt@ucdavis.edu")
        setTimeOut(david, milliSeconds = 10^6)
        setAnnotationCategories(david, categories = categories)
        addList(david, inputIds = genes, idType = "ENTREZ_GENE_ID", listName = "genes", listType = "Gene")[[1]]
        addList(david, inputIds = background, idType = "ENTREZ_GENE_ID", listName = "background", listType = "Background")[[1]]
        message("[getDAVID] Downloading enrichment results")
        results <- getFunctionalAnnotationChart(david) %>% subset(Benjamini < benjamini, select = -c(Bonferroni, FDR))
        results <- results[,c("Category", "Term", "Count", "Genes", "X.", "List.Total", "Pop.Hits", "Pop.Total", 
                              "Fold.Enrichment", "PValue", "Benjamini")]
        colnames(results) <- c("Category", "Term", "GeneListHits", "GeneIDs", "Percent", "GeneListTotal", "BackgroundHits",
                               "BackgroundTotal", "FoldEnrichment", "Pvalue", "BenjaminiPvalue")
        results$log10Benjamini <- -log10(results$BenjaminiPvalue)
        results <- results[order(results$BenjaminiPvalue, results$Pvalue),]
        message("[getDAVID] Complete! Writing file to ", file, "\n")
        write.table(results, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
        return(results)
}

prepGREAT <- function(DMRs, Background, fileName, writeBack = FALSE, backName, 
                      chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")){
        chain <- import.chain("Tables/hg38ToHg19.over.chain")
        
        message("[prepGREAT] Lifting over DMRs to hg19 and removing duplicates")
        seqlevelsStyle(DMRs) <- "UCSC"
        DMRs$RegionID <- paste("Region", 1:length(DMRs), sep = "_")
        DMRs_hg19 <- suppressWarnings(unlist(liftOver(DMRs, chain)))
        dups <- DMRs_hg19$RegionID[duplicated(DMRs_hg19$RegionID)] %>% unique
        DMRs_hg19_dups <- DMRs_hg19[DMRs_hg19$RegionID %in% dups]
        DMRs_hg19_unique <- DMRs_hg19[!DMRs_hg19$RegionID %in% dups]
        DMRs_hg19_dups_fix <- GRanges(NULL)
        for(i in 1:length(dups)){
                temp <- DMRs_hg19_dups[DMRs_hg19_dups$RegionID == dups[i]]
                temp <- temp[width(temp) == max(width(temp))] # Take largest duplicate
                DMRs_hg19_dups_fix <- union(DMRs_hg19_dups_fix, temp)
        }
        DMRs_hg19 <- union(DMRs_hg19_unique, DMRs_hg19_dups_fix)
        if(!isDisjoint(DMRs_hg19)){DMRs_hg19 <- disjoin(DMRs_hg19)} 
        
        message("[prepGREAT] Lifting over background to hg19 and removing duplicates")
        seqlevelsStyle(Background) <- "UCSC"
        Background$RegionID <- paste("Region", 1:length(Background), sep = "_")
        Background_hg19 <- suppressWarnings(unlist(liftOver(Background, chain)))
        dups <- Background_hg19$RegionID[duplicated(Background_hg19$RegionID)] %>% unique
        Background_hg19_dups <- Background_hg19[Background_hg19$RegionID %in% dups]
        Background_hg19_unique <- Background_hg19[!Background_hg19$RegionID %in% dups]
        Background_hg19_dups_fix <- GRanges(NULL)
        for(i in 1:length(dups)){
                temp <- Background_hg19_dups[Background_hg19_dups$RegionID == dups[i]]
                temp <- temp[width(temp) == max(width(temp))] # Take largest duplicate
                Background_hg19_dups_fix <- union(Background_hg19_dups_fix, temp)
        }
        Background_hg19 <- union(Background_hg19_unique, Background_hg19_dups_fix)
        if(!isDisjoint(Background_hg19)){Background_hg19 <- disjoin(Background_hg19)} 
        
        DMRs_hg19 <- redefineUserSets(GRangesList(DMRs_hg19), Background_hg19)
        message("[prepGREAT] ", table(overlapsAny(DMRs_hg19[[1]], Background_hg19, type = "equal")), 
                " DMRs in background out of ", length(overlapsAny(DMRs_hg19[[1]])), " total DMRs")
        
        message("[prepGREAT] Preparing DMR table")
        DMRs_hg19 <- as.data.frame(DMRs_hg19[[1]])[,c("seqnames", "start", "end")]
        colnames(DMRs_hg19) <- c("chr", "start", "end")
        DMRs_hg19$chr <- as.character(DMRs_hg19$chr)
        DMRs_hg19 <- DMRs_hg19[order(DMRs_hg19$chr, DMRs_hg19$start),]
        DMRs_hg19 <- unique(subset(DMRs_hg19, chr %in% chroms))
        write.table(DMRs_hg19, fileName, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
        message("[prepGREAT] Writing zipped DMR file to ", fileName)
        gzip(fileName, overwrite = TRUE)
        
        if(writeBack){
                message("[prepGREAT] Preparing background table")
                Background_hg19 <- as.data.frame(Background_hg19)[,c("seqnames", "start", "end")]
                colnames(Background_hg19) <- c("chr", "start", "end")
                Background_hg19$chr <- as.character(Background_hg19$chr)
                Background_hg19 <- Background_hg19[order(Background_hg19$chr, Background_hg19$start),]
                Background_hg19 <- unique(subset(Background_hg19, chr %in% chroms))
                if(nrow(Background_hg19) > 1000000){
                        message("[prepGREAT] Warning: Background has ", nrow(Background_hg19), 
                                " regions. Need to reduce background to < 1M regions")
                }
                write.table(Background_hg19, backName, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
                message("[prepGREAT] Writing zipped background file to ", backName)
                gzip(backName, overwrite = TRUE)
        }
        message("[prepGREAT] Complete!")
}

getGREATenrichments <- function(BEDfile_DMR, BEDfile_Back, species = "hg19", rule = "basalPlusExt", minGenes = 10, maxGenes = 2000,
                                exclude = c("MGI Expression: Detected", "MSigDB Oncogenic Signatures", 
                                            "MSigDB Cancer Neighborhood", "MSigDB Perturbation", 
                                            "MSigDB Predicted Promoter Motifs")){
        # Gets enrichments with FDR q < 0.05 from a DMR and Background bed.gz file, Fails with "MGI Expression: Detected"
        message("[getGREATenrichments] Submitting Regions to GREAT")
        job <- submitGreatJob(gr = BEDfile_DMR, bg = BEDfile_Back, species = species, includeCuratedRegDoms = FALSE, 
                              rule = rule, request_interval = 30, version = "3.0.0")
        message("[getGREATenrichments] Getting Enrichment Tables")
        allOntologies <- availableOntologies(job = job) %>% as.character
        enrichTables <- getEnrichmentTables(job = job, ontology = allOntologies[!allOntologies %in% exclude])
        message("[getGREATenrichments] Formatting Results")
        results <- NULL
        for(i in 1:length(enrichTables)){
                temp <- enrichTables[[i]]
                temp$Ontology <- rep(names(enrichTables)[i], nrow(temp))
                results <- rbind(results, temp)
        }
        results <- subset(results, Total_Genes_Annotated >= minGenes & Total_Genes_Annotated <= maxGenes)
        results$qvalue <- p.adjust(results$Hyper_Raw_PValue, method = "fdr")
        results <- subset(results, qvalue < 0.05)
        results$log_qvalue <- -log10(results$qvalue)
        results <- results[order(-results$log_qvalue),c("ID", "name", "Ontology", "Hyper_Foreground_Region_Hits",
                                                        "Hyper_Expected", "Hyper_Region_Set_Coverage", 
                                                        "Hyper_Term_Region_Coverage", "Hyper_Foreground_Gene_Hits",
                                                        "Hyper_Background_Gene_Hits", "Total_Genes_Annotated", 
                                                        "Hyper_Fold_Enrichment", "Hyper_Raw_PValue", "qvalue", "log_qvalue")]
        message("[getGREATenrichments] Complete!")
        return(results)
}

plotGREAT <- function(greatCombined, file, axis.text.y.size = 7.5, axis.text.y.width = 50, legend.position = c(1.32, 0.87),
                      width = 7, height = 5, wrap = FALSE){
        gg <- ggplot()
        gg <- gg +
                geom_tile(data = greatCombined, aes(y = name, x = Direction, fill = log_qvalue)) +
                scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                                     breaks = pretty_breaks(n = 3)) +
                scale_x_discrete(expand = c(0,0), drop = FALSE) +
                theme_bw(base_size = 24) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), 
                      plot.margin = unit(c(0.5, 5, 0.5, 0.5), "lines"), axis.ticks = element_line(size = 0.8), 
                      axis.text.x = element_text(size = 11, color = "black"), 
                      axis.text.y = element_text(size = axis.text.y.size, color = "black"),
                      axis.title = element_blank(), legend.key = element_blank(),  
                      legend.position = legend.position, legend.background = element_blank(), 
                      legend.key.size = unit(1, "lines"), legend.title = element_text(size = 11), 
                      legend.text = element_text(size = 11), panel.background = element_rect(fill = "black"))
        replacements <- c("Dna" = "DNA", "Of" = "of", "The" = "the", "Ensg" = "ENSG", "Ng" = "ng",
                          "Ml" = "ml", "Lps" = "LPS", "Pbmc" = "PBMC", "Tlr4" = "TLR4", 
                          "Trem1" = "TREM1", "With" = "with", " *\\(.*?\\)" = "", " *\\[.*?\\]" = "",
                          "\\." = "", "Genes " = "", "In " = "in ", "Rna" = "RNA", "Mrna" = "mRNA",
                          "Peripheral Blood Mononuclear Cells" = "PBMCs", "And" = "and", "Cd" = "CD",
                          "Nf" = "NF", "Kappab" = "kappaB", "Gtp" = "GTP", "To" = "to", "At " = "at ",
                          "Igg" = "IgG", "Cxcr" = "CXCR", "Ifn" = "IFN", "Dc" = "DC", "Tiv" = "TIV",
                          "Vaccinee" = "Vaccine", "Il4" = "IL4", "Ii" = "II", "NFat" = "NFAT", "Gnrh" = "GnRH",
                          "Er" = "ER", "Mtor" = "mTOR", "CDk" = "CDK", "Tgf" = "TGF", "Arf" = "ARF", "Mirna" = "miRNA",
                          "Il1" = "IL1", "ERythropoiesis" = "Erythropoiesis", "Hiv" = "HIV", "Hdac" = "HDAC", 
                          "Gaba" = "GABA", "Cns" = "CNS", "Hdl" = "HDL", "Microrna" = "MicroRNA", "Mir" = "MIR", 
                          "Cagtatt" = "CAGTATT", "Udp" = "UDP", "Hs-Gag" = "HS-GAG", "IIi" = "III", "Ch-Nh2" = "CH-NH2",
                          "androgen" = "Androgen", "ERrors" = "Errors", "Transcription Factor" = "TF")
        if(wrap){
                gg <- gg +
                        scale_y_discrete(expand = c(0,0), drop = FALSE, labels = function(x){
                                str_to_title(x) %>% str_replace_all(replacements) %>% 
                                        str_wrap(width = axis.text.y.width) %>% 
                                        str_trunc(width = axis.text.y.width * 2 - 8, side = "right")
                        })
        }
        else {
                gg <- gg +
                        scale_y_discrete(expand = c(0,0), drop = FALSE, labels = function(x){
                                str_to_title(x) %>% str_replace_all(replacements) %>% 
                                        str_trunc(width = axis.text.y.width * 2 - 8, side = "right")
                        })
        }
        ggsave(file, plot = gg, dpi = 600, width = width, height = height, units = "in")
}

prepLOLAhistone <- function(histone, index, regions, file){
        message("[prepLOLAhistone] Preparing histone enrichment results for ", regions)
        histone <- subset(histone, userSet == regions & antibody %in% c("H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"))
        histone$EID <- strsplit(as.character(histone$filename), "-") %>% sapply(function(x){x[1]})
        match <- match(histone$EID, index$EID)
        histone$type <- index$type[match]
        histone$order <- index$order[match]
        histone$color <- index$color[match]
        histone$cellType <- index$cellType[match]
        histone$tissue <- index$tissue[match]
        
        if(!is.na(table(histone$support < 5)["TRUE"]) & table(histone$support < 5)["TRUE"] > 0){
                message("[prepLOLAhistone] Editing overlaps with support < 5")
                histone$oddsRatio[histone$support < 5] <- 0
                histone$pValueLog[histone$support < 5] <- 0
                histone$qValue[histone$support < 5] <- 1
        } 
        histone$pct_DMRs <- histone$support * 100 / (histone$support[1] + histone$c[1])
        
        if(!is.na(table(is.infinite(histone$oddsRatio))["TRUE"]) & table(is.infinite(histone$oddsRatio))["TRUE"] > 0){
                message("[prepLOLAhistone] Replacing infinite odds ratios with max odds ratio for that antibody")
                histone$oddsRatio[is.infinite(histone$oddsRatio)] <- NA
                replace <- unique(as.character(histone$antibody[is.na(histone$oddsRatio)]))
                for(i in 1:length(replace)){
                        histone$oddsRatio[histone$antibody == replace[i] & is.na(histone$oddsRatio)] <- 
                                max(histone$oddsRatio[histone$antibody == replace[i]], na.rm = TRUE)
                }
        }
        
        if(!is.na(table(is.infinite(histone$pValueLog))["TRUE"]) & table(is.infinite(histone$pValueLog))["TRUE"] > 0){
                message("[prepLOLAhistone] Replacing infinite log(p-values) with max")
                histone$pValueLog[is.infinite(histone$pValueLog)] <- NA
                histone$pValueLog[is.na(histone$pValueLog)] <- max(histone$pValueLog, na.rm = TRUE)
        }
        
        histone$pValue <- 10^(-histone$pValueLog)
        histone$qValueLog <- -log10(histone$qValue)
        if(!is.na(table(is.infinite(histone$qValueLog))["TRUE"]) & table(is.infinite(histone$qValueLog))["TRUE"] > 0){
                message("[prepLOLAhistone] Replacing infinite log(q-values) with max")
                histone$qValueLog[is.infinite(histone$qValueLog)] <- NA
                histone$qValueLog[is.na(histone$qValueLog)] <- max(histone$qValueLog, na.rm = TRUE)
        }
        
        histone <- histone[,c("userSet", "dbSet", "antibody", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", 
                              "oddsRatio", "support", "pct_DMRs", "rnkPV", "rnkOR", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", 
                              "size", "filename", "order", "color")]
        histone$order <- factor(histone$order, levels = sort(unique(histone$order), decreasing = TRUE), ordered = TRUE)
        histone$antibody <- factor(histone$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"), 
                                   ordered = TRUE)
        histone <- histone[order(histone$order, histone$antibody), ]
        histone$tissue <- factor(histone$tissue, levels = rev(as.character(unique(histone$tissue))), ordered = TRUE)
        message("[prepLOLAhistone] Complete! Writing file")
        write.csv(x = histone, file = file, quote = FALSE, row.names = FALSE)
        return(histone)
}

plotLOLAhistone <- function(histone, title, type = c("oddsRatio", "qValueLog", "legend"), hm.max, file, axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(), labels = NULL, width = 5.5, height = 7, legend.position = c(1.21, 0.855), 
                            facet = NULL){
        if(!type %in% c("oddsRatio", "qValueLog", "legend")){
                message("[plotLOLAhistone] type must be oddsRatio, qValueLog or legend")
        }
        else {
                message("[plotLOLAhistone] Plotting histone enrichment ", type)
                if(type == "oddsRatio"){
                        if(!is.null(facet)){
                                gg <- ggplot(data = histone)
                                gg +
                                        geom_tile(aes(x = antibody, y = order, fill = oddsRatio)) +
                                        facet_grid(cols = facet) +
                                        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(),  legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        } else {
                                gg <- ggplot(data = histone)
                                gg +
                                        geom_tile(aes(x = antibody, y = order, fill = oddsRatio)) +
                                        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(),  legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        }
                }
                if(type == "qValueLog"){
                        if(!is.null(facet)){
                                gg <- ggplot(data = histone)
                                gg +
                                        geom_tile(aes(x = antibody, y = order, fill = qValueLog)) +
                                        facet_grid(cols = facet) +
                                        scale_fill_gradientn("-log(q-value)", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(), legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        } else {
                                gg <- ggplot(data = histone)
                                gg +
                                        geom_tile(aes(x = antibody, y = order, fill = qValueLog)) +
                                        scale_fill_gradientn("-log(q-value)", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(), legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        }
                }
                if(type == "legend"){
                        gg <- ggplot(data = histone)
                        gg +
                                geom_tile(aes(x = 16, y = order, fill = tissue)) +
                                scale_fill_manual(name = "Tissue", values = rev(as.character(unique(histone$color)))) +
                                theme_bw(base_size = 24) +
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                      panel.border = element_rect(color = "black", size = 1.25), legend.key = element_blank(), 
                                      legend.position = c(4.7, 0.23), legend.background = element_blank(), 
                                      legend.text = element_text(size = 15, color = "Black"), plot.margin = unit(c(2.5, 22, 6.3, 1), "lines"), 
                                      axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
                                      legend.title = element_text(size = 18)) +
                                scale_x_discrete(expand = c(0, 0))
                        ggsave(file, dpi = 600, width = 5, height = 7, units = "in")
                }
                message("[plotLOLAhistone] Complete! Writing file")
        }
}

plotLOLAhistoneBox <- function(histone, file, facet = NULL, breaks = c("Autosomes", "ChrX"), 
                               values = c("Autosomes" = "#3366CC", "ChrX" = "#FF3366"), plot.margin = c(1.5, 0.5, 0.5, 1.15),
                               legend.position = c(0.82, 1.12), ylim = NULL){
        gg <- ggplot()
        gg +
                geom_abline(slope = 0, intercept = 1, lwd = 1.25, lty = 2) +
                stat_summary(data = histone, aes(x = antibody, y = oddsRatio, color = chroms), fun.data = "mean_cl_boot",
                             geom = "crossbar", size = 0.7, position = position_dodge2(width = 0.5, padding = 0.15), fill = "white") +
                theme_bw(base_size = 24) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
                      legend.key = element_blank(), legend.position = legend.position, legend.text = element_text(size = 20),
                      legend.background = element_blank(), legend.title = element_blank(), legend.spacing.x = unit(0.5, "lines"),
                      plot.margin = unit(plot.margin, "lines"), strip.text = element_text(size = 20), 
                      axis.text.x = element_text(size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5),
                      axis.text.y = element_text(size = 18, color = "black"), axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 20, color = "black"),
                      strip.background = element_blank()) +
                scale_color_manual(breaks = breaks, values = values) +
                scale_y_continuous(breaks = pretty_breaks(6)) +
                coord_cartesian(ylim = ylim) +
                ylab("Odds Ratio") +
                facet_grid(cols = facet, scales = "fixed")
        ggsave(file, dpi = 600, width = 9, height = 7, units = "in")
}

prepLOLAchromHMM <- function(chromHMM, index, regions, file){
        message("[prepLOLAchromHMM] Preparing chromHMM enrichment results for ", regions)
        chromHMM <- subset(chromHMM, userSet == regions)
        chromHMM$EID <- strsplit(as.character(chromHMM$filename), "_") %>% sapply(function(x){x[1]})
        match <- match(chromHMM$EID, index$EID)
        chromHMM$type <- index$type[match]
        chromHMM$order <- index$order[match]
        chromHMM$color <- index$color[match]
        chromHMM$cellType <- index$cellType[match]
        chromHMM$tissue <- index$tissue[match]
        
        if(!is.na(table(chromHMM$support < 5)["TRUE"]) & table(chromHMM$support < 5)["TRUE"] > 0){
                message("[prepLOLAchromHMM] Editing overlaps with support < 5")
                chromHMM$oddsRatio[chromHMM$support < 5] <- 0
                chromHMM$pValueLog[chromHMM$support < 5] <- 0
                chromHMM$qValue[chromHMM$support < 5] <- 1
        } 
        chromHMM$pct_DMRs <- chromHMM$support * 100 / (chromHMM$support[1] + chromHMM$c[1])
        
        if(!is.na(table(is.infinite(chromHMM$oddsRatio))["TRUE"]) & table(is.infinite(chromHMM$oddsRatio))["TRUE"] > 0){
                message("[prepLOLAchromHMM] Replacing infinite odds ratios with max odds ratio for that antibody")
                chromHMM$oddsRatio[is.infinite(chromHMM$oddsRatio)] <- NA
                replace <- unique(as.character(chromHMM$antibody[is.na(chromHMM$oddsRatio)]))
                for(i in 1:length(replace)){
                        chromHMM$oddsRatio[chromHMM$antibody == replace[i] & is.na(chromHMM$oddsRatio)] <- 
                                max(chromHMM$oddsRatio[chromHMM$antibody == replace[i]], na.rm = TRUE)
                }
        }
        
        if(!is.na(table(is.infinite(chromHMM$pValueLog))["TRUE"]) & table(is.infinite(chromHMM$pValueLog))["TRUE"] > 0){
                message("[prepLOLAchromHMM] Replacing infinite log(p-values) with max")
                chromHMM$pValueLog[is.infinite(chromHMM$pValueLog)] <- NA
                chromHMM$pValueLog[is.na(chromHMM$pValueLog)] <- max(chromHMM$pValueLog, na.rm = TRUE)
        }
        
        chromHMM$pValue <- 10^(-chromHMM$pValueLog)
        chromHMM$qValueLog <- -log10(chromHMM$qValue)
        if(!is.na(table(is.infinite(chromHMM$qValueLog))["TRUE"]) & table(is.infinite(chromHMM$qValueLog))["TRUE"] > 0){
                message("[prepLOLAchromHMM] Replacing infinite log(q-values) with max")
                chromHMM$qValueLog[is.infinite(chromHMM$qValueLog)] <- NA
                chromHMM$qValueLog[is.na(chromHMM$qValueLog)] <- max(chromHMM$qValueLog, na.rm = TRUE)
        }
        
        chromHMM <- chromHMM[,c("userSet", "dbSet", "antibody", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", 
                                "oddsRatio", "support", "pct_DMRs", "rnkPV", "rnkOR", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", 
                                "size", "filename", "order", "color")]
        chromHMM$order <- factor(chromHMM$order, levels = sort(unique(chromHMM$order), decreasing = TRUE), ordered = TRUE)
        colnames(chromHMM)[colnames(chromHMM) == "antibody"] <- "chromState"
        chromHMM$chromState <- sapply(strsplit(as.character(chromHMM$chromState), split='_', fixed=TRUE), function(x) (x[2]))
        chromHMM$chromState <- factor(chromHMM$chromState, levels = c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", 
                                                                      "ZnfRpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", 
                                                                      "ReprPCwk", "Quies"), ordered = TRUE)
        chromHMM <- chromHMM[order(chromHMM$order, chromHMM$chromState),]
        chromHMM$tissue <- factor(chromHMM$tissue, levels = rev(as.character(unique(chromHMM$tissue))), ordered = TRUE)
        message("[prepLOLAchromHMM] Complete! Writing file")
        write.csv(x = chromHMM, file = file, quote = FALSE, row.names = FALSE)
        return(chromHMM)
}

plotLOLAchromHMM <- function(chromHMM, title, type = c("oddsRatio", "qValueLog", "legend"), hm.max, file, axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(), labels = NULL, width = 5.5, height = 7, legend.position = c(1.21, 0.855), 
                            facet = NULL){
        if(!type %in% c("oddsRatio", "qValueLog", "legend")){
                message("[plotLOLAchromHMM] type must be oddsRatio, qValueLog or legend")
        }
        else {
                message("[plotLOLAchromHMM] Plotting chromHMM enrichment ", type)
                if(type == "oddsRatio"){
                        if(!is.null(facet)){
                                gg <- ggplot(data = chromHMM)
                                gg +
                                        geom_tile(aes(x = chromState, y = order, fill = oddsRatio)) +
                                        facet_grid(cols = facet) +
                                        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(),  legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        } else {
                                gg <- ggplot(data = chromHMM)
                                gg +
                                        geom_tile(aes(x = chromState, y = order, fill = oddsRatio)) +
                                        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(),  legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        }
                }
                if(type == "qValueLog"){
                        if(!is.null(facet)){
                                gg <- ggplot(data = chromHMM)
                                gg +
                                        geom_tile(aes(x = chromState, y = order, fill = qValueLog)) +
                                        facet_grid(cols = facet) +
                                        scale_fill_gradientn("-log(q-value)", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(), legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        } else {
                                gg <- ggplot(data = chromHMM)
                                gg +
                                        geom_tile(aes(x = chromState, y = order, fill = qValueLog)) +
                                        scale_fill_gradientn("-log(q-value)", colors = c("black", "#FF0000"), values = c(0, 1), 
                                                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
                                        theme_bw(base_size = 24) +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                              panel.border = element_rect(color = "black", size = 1.25), 
                                              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
                                              axis.ticks.y = axis.ticks.y, legend.key = element_blank(), legend.position = legend.position, 
                                              legend.background = element_blank(), legend.title = element_text(size = 18), 
                                              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = axis.text.y, 
                                              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
                                              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                                              strip.background = element_blank(), strip.text = element_text(size = 18)) +
                                        scale_x_discrete(expand = c(0, 0)) + 
                                        scale_y_discrete(labels = labels) +
                                        ggtitle(title)
                                ggsave(file, dpi = 600, width = width, height = height, units = "in")
                        }
                }
                if(type == "legend"){
                        gg <- ggplot(data = chromHMM)
                        gg +
                                geom_tile(aes(x = 16, y = order, fill = tissue)) +
                                scale_fill_manual(name = "Tissue", values = rev(as.character(unique(chromHMM$color)))) +
                                theme_bw(base_size = 24) +
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                      panel.border = element_rect(color = "black", size = 1.25), legend.key = element_blank(), 
                                      legend.position = c(4.7, 0.23), legend.background = element_blank(), 
                                      legend.text = element_text(size = 15, color = "Black"), plot.margin = unit(c(2.5, 22, 6.3, 1), "lines"), 
                                      axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
                                      legend.title = element_text(size = 18)) +
                                scale_x_discrete(expand = c(0, 0))
                        ggsave(file, dpi = 600, width = 5, height = 7, units = "in")
                }
                message("[plotLOLAchromHMM] Complete! Writing file")
        }
}

plotLOLAchromHMMBox <- function(chromHMM, file, facet = NULL, breaks = c("Autosomes", "ChrX"), 
                                values = c("Autosomes" = "#3366CC", "ChrX" = "#FF3366"), plot.margin = c(1.5, 0.5, 0.75, 0.5),
                                legend.position = c(0.82, 1.12), ylim = NULL){
        gg <- ggplot()
        gg +
                geom_abline(slope = 0, intercept = 1, lwd = 1.25, lty = 2) +
                stat_summary(data = chromHMM, aes(x = chromState, y = oddsRatio, color = chroms), fun.data = "mean_cl_boot",
                             geom = "crossbar", size = 0.7, position = position_dodge2(width = 0.5, padding = 0.15), fill = "white") +
                theme_bw(base_size = 24) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.25),
                      axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
                      legend.key = element_blank(), legend.position = legend.position, legend.text = element_text(size = 20),
                      legend.background = element_blank(), legend.title = element_blank(), legend.spacing.x = unit(0.5, "lines"),
                      plot.margin = unit(plot.margin, "lines"), strip.text = element_text(size = 20), 
                      axis.text.x = element_text(size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5),
                      axis.text.y = element_text(size = 18, color = "black"), axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 20, color = "black"),
                      strip.background = element_blank()) +
                scale_color_manual(breaks = breaks, values = values) +
                scale_y_continuous(breaks = pretty_breaks(6)) +
                coord_cartesian(ylim = ylim) +
                ylab("Odds Ratio") +
                facet_grid(cols = facet, scales = "fixed")
        ggsave(file, dpi = 600, width = 9, height = 7, units = "in")
}

loadHOMER <- function(file){
        homer <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% unique()
        homer <- homer[order(homer$Motif.Name, homer$Consensus),]
        Motif.Name <- strsplit(homer$Motif.Name, split = "/", fixed = TRUE) %>%
                sapply(function(x) strsplit(x[1], split = "(", fixed = TRUE)[[1]])
        homer2 <- data.frame(stringsAsFactors = FALSE,
                             Transcription_Factor = sapply(Motif.Name, function(x) x[1]) %>% str_to_upper() %>% 
                                     str_replace_all(pattern = fixed(c("UNKNOWN-ESC-ELEMENT" = "Unknown_ESC_Element", 
                                                                       "DISTAL" = "Distal", "FUSION" = "Fusion", 
                                                                       "UNKNOWN" = "Unknown", "SHORT" = "Short", 
                                                                       "HALFSITE" = "Half_Site", "|" = ":",
                                                                       "MOUSE_RECOMBINATION_HOTSPOT" = "Mouse_Recombination_Hotspot", 
                                                                       "SATELLITEELEMENT" = "Satellite_Element", "?" = "", 
                                                                       "BRACHYURY" = "Brachyury"))),
                             Family = sapply(Motif.Name, function(x) x[2]) %>% 
                                     str_replace_all(pattern = fixed(c(")" = "", "?," = "", "?" = "", "," = "_"))),
                             Consensus = homer$Consensus,
                             Target_with_Motif = homer[,grepl("X..of.Target.Sequences.with.Motif.of." , x = colnames(homer), 
                                                              fixed = TRUE)],
                             Target_Sequences = colnames(homer)[grepl("X..of.Target.Sequences.with.Motif.of." , 
                                                                      x = colnames(homer), fixed = TRUE)] %>% 
                                     strsplit(split = ".", fixed = TRUE) %>% .[[1]] %>% .[length(.)] %>% as.numeric(),
                             Background_with_Motif = homer[,grepl("X..of.Background.Sequences.with.Motif.of." , 
                                                                  x = colnames(homer), fixed = TRUE)],
                             Background_Sequences = colnames(homer)[grepl("X..of.Background.Sequences.with.Motif.of." , 
                                                                          x = colnames(homer), fixed = TRUE)] %>% 
                                     strsplit(split = ".", fixed = TRUE) %>% .[[1]] %>% .[length(.)] %>% as.numeric(),
                             pvalue = homer$P.value)
        dups <- homer2$Transcription_Factor[duplicated(homer2$Transcription_Factor)] %>% unique()
        for(i in 1:length(dups)){
                homer2$Transcription_Factor[homer2$Transcription_Factor == dups[i]] <- 
                        paste(homer2$Transcription_Factor[homer2$Transcription_Factor == dups[i]],
                              1:length(homer2$Transcription_Factor[homer2$Transcription_Factor == dups[i]]),
                              sep = "_")
        }
        homer2$Family[homer2$Family == ""] <- "Unknown"
        homer2$Percent_Target_with_Motif <- homer2$Target_with_Motif * 100 / homer2$Target_Sequences
        homer2$Percent_Background_with_Motif <- homer2$Background_with_Motif * 100 / homer2$Background_Sequences
        homer2$Fold_Enrichment <- homer2$Percent_Target_with_Motif / homer2$Percent_Background_with_Motif
        homer2$log_pvalue <- -log10(homer2$pvalue)
        homer2$qvalue <- p.adjust(homer2$pvalue, method = "fdr")
        homer2$log_qvalue <- -log10(homer2$qvalue)
        homer2 <- homer2[order(homer2$log_qvalue, homer2$Fold_Enrichment, decreasing = TRUE),
                         c("Transcription_Factor", "Family", "Consensus", "Target_with_Motif", "Target_Sequences", 
                           "Percent_Target_with_Motif", "Background_with_Motif", "Background_Sequences",
                           "Percent_Background_with_Motif", "Fold_Enrichment", "pvalue", "log_pvalue", "qvalue", 
                           "log_qvalue")]
        return(homer2)
}

# Replication Functions ---------------------------------------------------

DMRoverlapVenn <- function(Peaks, NameOfPeaks, file, totalTest = NULL, rotation.degree = 0, cat.pos = c(0, 0), 
                           cat.dist = c(0.03, 0.03), cat.cex = 3, fill = c("lightblue", "lightpink"), ext.text = TRUE, 
                           ext.dist = -0.2, ext.percent = c(0.01, 0.01, 0.01), margin = 0.04, cex = 2.5){
        pdf(file = file, width = 10, height = 8, onefile = FALSE)
        if(!is.null(totalTest)){
                venn <- suppressMessages(makeVennDiagram(Peaks = Peaks, NameOfPeaks = NameOfPeaks, totalTest = totalTest, 
                                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                                         rotation.degree = rotation.degree, margin = margin, cat.cex = cat.cex, 
                                                         cex = cex, fill = fill, cat.pos = cat.pos, cat.dist = cat.dist, 
                                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = ext.dist, 
                                                         ext.length = 0.85, ext.text = ext.text, ext.percent = ext.percent, 
                                                         lwd = 4, ext.line.lwd = 3))
                message("[DMRoverlapVenn] Hypergeometric test p = ", signif(venn$p.value[3], digits = 6), "\n")
        } else {
                venn <- suppressMessages(makeVennDiagram(Peaks = Peaks, NameOfPeaks = NameOfPeaks, maxgap = -1, minoverlap = 1, 
                                                         by = "region", connectedPeaks = "min", rotation.degree = rotation.degree, 
                                                         margin = margin, cat.cex = cat.cex, cex = cex, fill = fill,
                                                         cat.pos = cat.pos, cat.dist = cat.dist, fontfamily = "sans",
                                                         cat.fontfamily = "sans", ext.dist = ext.dist, ext.length = 0.85, 
                                                         ext.text = ext.text, ext.percent = ext.percent, lwd = 4,
                                                         ext.line.lwd = 3))
        }
        dev.off()
}

geneOverlapVenn <- function(x, file, cat.pos = c(150, 180), cat.dist = c(0.04, 0.03), rotation.degree = 180,
                            ext.dist = -0.1, margin = 0.05, cat.cex = 2.5, cex = 2.75, fill = c("lightblue", "lightpink"), reverse = FALSE){
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # Suppress log file
        venn.diagram(x = x, file = file, height = 8, width = 10, imagetype = "png", units = "in", fontfamily = "sans", 
                     cat.fontfamily = "sans", fill = fill, cex = cex, lwd = 4, cat.cex = cat.cex, 
                     cat.pos = cat.pos, cat.dist = cat.dist, rotation.degree = rotation.degree, margin = margin, 
                     ext.text = TRUE, ext.dist = ext.dist, ext.length = 0.85, ext.line.lwd = 3, 
                     ext.percent = 0.01, reverse = reverse)
}

plotGOMheatmap <- function(data, type = c("log_OddsRatio", "log_qValue"), file, expand = c(0.04, 0), limits = NULL,
                           intersect.size = 5, sig.size = 9, axis.text.size = 14, legend.text.size = 14, legend.title.size = 16,
                           plot.margin = c(1, 8, 1, 1), legend.position = c(1.09, 0.87), width = 10, height = 8){
        if(type == "log_OddsRatio"){
                gg <- ggplot(data = data)
                gg +
                        geom_tile(aes(x = ListA, y = ListB, fill = log10(OddsRatio))) +
                        geom_text(aes(x = ListA, y = ListB, label = Intersection), color = "white", size = intersect.size, nudge_y = -0.15) +
                        geom_text(aes(x = ListA, y = ListB, label = "*", alpha = Significant), color = "white", size = sig.size, nudge_y = 0.05) +
                        scale_fill_gradientn("log(OR)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "black", limits = limits) +
                        scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0), guide = "none") +
                        theme_bw(base_size = 24) +
                        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
                              axis.ticks = element_line(size = 1), legend.key = element_blank(), 
                              legend.text = element_text(size = legend.text.size),
                              panel.grid.minor = element_blank(), legend.position = legend.position, 
                              legend.background = element_blank(), panel.background = element_rect(fill = "#333333"),
                              plot.margin = unit(plot.margin, "lines"), 
                              axis.text.x = element_text(size = axis.text.size, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
                              axis.text.y = element_text(size = axis.text.size, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
                              axis.title = element_blank(), legend.title = element_text(size = legend.title.size)) +
                        scale_x_discrete(expand = expand, drop = FALSE) +
                        scale_y_discrete(expand = expand, drop = FALSE)
                ggsave(file = file, dpi = 600, width = width, height = height, units = "in")
        }
        else {
                gg <- ggplot(data = data)
                gg +
                        geom_tile(aes(x = ListA, y = ListB, fill = log_qValue)) +
                        geom_text(aes(x = ListA, y = ListB, label = Intersection), color = "white", size = intersect.size, nudge_y = -0.15) +
                        geom_text(aes(x = ListA, y = ListB, label = "*", alpha = Significant), color = "white", size = sig.size, nudge_y = 0.05) +
                        scale_fill_gradientn("log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "black", limits = limits) +
                        scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0), guide = "none") +
                        theme_bw(base_size = 24) +
                        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
                              axis.ticks = element_line(size = 1), legend.key = element_blank(), legend.text = element_text(size = legend.text.size),
                              panel.grid.minor = element_blank(), legend.position = legend.position, 
                              legend.background = element_blank(), panel.background = element_rect(fill = "#333333"),
                              plot.margin = unit(plot.margin, "lines"), 
                              axis.text.x = element_text(size = axis.text.size, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
                              axis.text.y = element_text(size = axis.text.size, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
                              axis.title = element_blank(), legend.title = element_text(size = legend.title.size)) +
                        scale_x_discrete(expand = expand, drop = FALSE) +
                        scale_y_discrete(expand = expand, drop = FALSE)
                ggsave(file = file, dpi = 600, width = width, height = height, units = "in")
        }
}

DMRpermTest <- function(A, B, genome, universe, Comparison, file, ntimes = 10000){
        message("[DMRpermTest] Performing permutation test of regions using regioneR.")
        pt <- permTest(A = A, B = B, genome = genome, ntimes = ntimes, universe = universe, 
                       evaluate.function = c(numOverlaps, meanDistance), randomize.function = resampleRegions, 
                       mc.set.seed = FALSE, force.parallel = TRUE)
        stats <- data.frame("Comparison" = Comparison, "Overlap_observed" = pt$Function1$observed, 
                            "Overlap_zscore" = pt$Function1$zscore, "Overlap_pvalue" = pt$Function1$pval, 
                            "Distance_observed" = pt$Function2$observed, "Distance_zscore" = pt$Function2$zscore, 
                            "Distance_pvalue" = pt$Function2$pval)
        message("[DMRpermTest] Complete! Writing plot and returning stats.")
        pdf(file = file, width = 10, height = 5)
        plot(x = pt, ncol = 2)
        dev.off()
        return(stats)
}

DMRrandomForest <- function(discDMRmeth, repDMRmeth, samples){
        # Prepare Data
        message("Preparing Discovery Methylation Data")
        discMeth <- discDMRmeth[,grepl("JLCM", colnames(discDMRmeth), fixed = TRUE)] %>% t %>% as.data.frame
        colnames(discMeth) <- discDMRmeth$DMRid
        if(table(is.na(discMeth))["TRUE"] > 0){
                message("Replacing NA values with DMR mean methylation")
                for(i in 1:ncol(discMeth)){
                        temp <- as.numeric(discMeth[,i])
                        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with mean meth of that DMR
                        discMeth[,i] <- temp
                }
        }
        discMeth$Diagnosis <- factor(samples$Diagnosis_Alg[match(rownames(discMeth), samples$Sequencing_ID)], 
                                     levels = c("TD", "ASD"))
        
        message("Preparing Replication Discovery Methylation Data")
        repMeth <- repDMRmeth[,grepl("JLCM", colnames(repDMRmeth), fixed = TRUE)] %>% t %>% as.data.frame
        colnames(repMeth) <- repDMRmeth$DMRid
        if(table(is.na(repMeth))["TRUE"] > 0){
                message("Replacing NA values with DMR mean methylation")
                for(i in 1:ncol(repMeth)){
                        temp <- as.numeric(repMeth[,i])
                        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with mean meth of that DMR
                        repMeth[,i] <- temp
                }
        }
        repMeth$Diagnosis <- factor(samples$Diagnosis_Alg[match(rownames(repMeth), samples$Sequencing_ID)], 
                                    levels = c("TD", "ASD"))
        
        # Initial Model
        message("\nMaking initial random forest model", appendLF = FALSE)
        set.seed(5)
        model1 <- randomForest(Diagnosis ~ ., data = discMeth, localImp = TRUE, ntree = 500)
        print(model1)
        
        message("Predicting discovery group diagnosis with initial model")
        discInitialPredict <- predict(model1, discMeth, type = "class")
        message("Classifications", appendLF = FALSE)
        print(table(discInitialPredict, discMeth$Diagnosis))  
        message("Accuracy: ", mean(discInitialPredict == discMeth$Diagnosis))
        
        message("\nPredicting replication group diagnosis with initial model")
        repInitialPredict <- predict(model1, repMeth, type = "class")
        message("Classifications", appendLF = FALSE)
        print(table(repInitialPredict,repMeth$Diagnosis))
        message("Accuracy: ", mean(repInitialPredict == repMeth$Diagnosis))
        
        # Optimize parameters
        message("\nOptimizing parameters")
        message("mtry =\t", appendLF = FALSE)
        def <- nrow(discDMRmeth) %>% sqrt %>% round
        mtryOpt <- NULL
        for(mtry in round(def/4):(def*4)){
                if(mtry %% 10 == 0){message(mtry, "\t", appendLF = FALSE)}
                set.seed(5)
                model <- randomForest(Diagnosis ~ ., data = discMeth, mtry = mtry, ntree = 500, localImp = TRUE)
                error <- tail(model$err.rate[,"OOB"], 1) * 100
                temp <- c(mtry, 500, error)
                mtryOpt <- rbind(mtryOpt, temp)
        }
        mtryOpt <- as.data.frame(mtryOpt)
        rownames(mtryOpt) <- 1:nrow(mtryOpt)
        colnames(mtryOpt) <- c("mtry", "ntree", "OOBerror")
        message("\nLowest OOB Error: ", min(mtryOpt$OOBerror))
        mtryOpt <- subset(mtryOpt, OOBerror == min(mtryOpt$OOBerror))
        print(mtryOpt)
        mtryOpt <- min(mtryOpt$mtry)
        message("Optimized mtry value: ", mtryOpt)
        
        message("ntree =\t", appendLF = FALSE)
        ntreeOpt <- NULL
        for(ntree in seq(100, 2500, 50)){
                if(ntree %% 500 == 0){message(ntree, "\t", appendLF = FALSE)}
                set.seed(5)
                model <- randomForest(Diagnosis ~ ., data = discMeth, mtry = mtryOpt, ntree = ntree, localImp = TRUE)
                error <- tail(model$err.rate[,"OOB"], 1) * 100
                temp <- c(mtryOpt, ntree, error)
                ntreeOpt <- rbind(ntreeOpt, temp)
        }
        ntreeOpt <- as.data.frame(ntreeOpt)
        rownames(ntreeOpt) <- 1:nrow(ntreeOpt)
        colnames(ntreeOpt) <- c("mtry", "ntree", "OOBerror")
        message("\nLowest OOB Error: ", min(ntreeOpt$OOBerror))
        ntreeOpt <- subset(ntreeOpt, OOBerror == min(ntreeOpt$OOBerror))
        print(ntreeOpt)
        ntreeOpt <- min(ntreeOpt$ntree)
        message("Optimized ntree value: ", ntreeOpt)
        
        # Optimized Model
        message("\nMaking optimized random forest model")
        set.seed(5)
        model2 <- randomForest(Diagnosis ~ ., data = discMeth, localImp = TRUE, mtry = mtryOpt, ntree = ntreeOpt)
        print(model2)
        
        message("Predicting discovery group diagnosis with optimized model")
        discOptimizedPredict <- predict(model2, discMeth, type = "class")
        message("Classifications", appendLF = FALSE)
        print(table(discOptimizedPredict, discMeth$Diagnosis))
        message("Accuracy: ", mean(discOptimizedPredict == discMeth$Diagnosis))
        
        message("\nPredicting replication group diagnosis with optimized model")
        repOptimizedPredict <- predict(model2, repMeth, type = "class")
        message("Classifications", appendLF = FALSE)
        print(table(repOptimizedPredict,repMeth$Diagnosis))
        message("Accuracy: ", mean(repOptimizedPredict == repMeth$Diagnosis))
        
        discPredict <- data.frame("InitialModelPredict" = discInitialPredict, "OptimizedModelPredict" = discOptimizedPredict, 
                                  "Diagnosis" = discMeth$Diagnosis)
        repPredict <- data.frame("InitialModelPredict" = repInitialPredict, "OptimizedModelPredict" = repOptimizedPredict,
                                 "Diagnosis" = repMeth$Diagnosis)
        
        # Results
        message("\n Returning model information")
        modelInfo <- list("discMeth" = discMeth, "repMeth" = repMeth, "InitialModel" = model1, "OptimizedModel" = model2,
                          "discPredict" = discPredict, "repPredict" = repPredict)
        
        message("Initial Model Summary")
        print(modelInfo$InitialModel)
        message("Classifications")
        print(table(modelInfo$repPredict$InitialModelPredict, modelInfo$repPredict$Diagnosis))
        message("Accuracy = ", mean(modelInfo$repPredict$InitialModelPredict == modelInfo$repPredict$Diagnosis))
        
        message("\nOptimized Model Summary")
        print(modelInfo$OptimizedModel)
        message("Classifications")
        print(table(modelInfo$repPredict$OptimizedModelPredict, modelInfo$repPredict$Diagnosis))
        message("Accuracy = ", mean(modelInfo$repPredict$OptimizedModelPredict == modelInfo$repPredict$Diagnosis))
        return(modelInfo)
}


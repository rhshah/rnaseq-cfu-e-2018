---
title: '&quot;Clustering and differential analysis of CFU-E cells with different conditions&quot;'
author: "Ronak Shah"
date: "July 27th 2017"
output:
  html_document:
    code_download: yes
    code_folding: hide
    highlight: pygments
    keep_md: yes
    theme: paper
    toc: yes
    toc_float: yes
    toc_depth: 5
  word_document:
    toc: yes
  pdf_document: default
---



## Table of Contents

- [Summary of Experiment:](#summary-of-experiment:)
- [Exploring data with DESeq2](#exploring-data-with-deseq2)
    - [Running DESeq2](#running-deseq2)
    - [Make counts table to explore that data](#make-counts-table-to-explore-that-data)
    - [Histograms of counts per gene](#histograms-of-counts-per-gene)
        - [matrix](#matrix)
        - [truncated matrix](#truncated-matrix)
        - [box-plotdensity-plots](#box-plotdensity-plots)
        - [scatter-plots](#scatter-plots)
    - [Eliminate undetected genes](#eliminate-undetected-genes)
    - [ Re-run DESeq2 & Normalization ](#-re-run-deseq2-&-normalization)
        - [Plot normalized data](#plot-normalized-data)
    - [Count variance is related to mean](#count-variance-is-related-to-mean)
    - [Estimated Dispersion for each gene](#estimated-dispersion-for-each-gene)
    - [PCA and heatmap on all samples using reguralized log transformation](#pca-and-heatmap-on-all-samples-using-reguralized-log-transformation)
        - [PCA using rlog](#pca-using-rlog)
        - [Heatmap of count matrix for first 50 genes](#heatmap-of-count-matrix-for-first-50-genes)
        - [Heatmap of sample distance based on rlog](#heatmap-of-sample-distance-based-on-rlog)
    - [Perform differential expression call](#perform-differential-expression-call)
        - [Set standard cutoff and Enhanced Volcano Plots](#set-standard-cutoff-and-enhanced-volcano-plots)
        - [Untreated (WT) vs EPO with Vehicle (EV) ](#untreated-(wt)-vs-epo-with-vehicle-(ev))
        - [Untreated (WT) vs EPO with HMGB1 (EH)](#untreated-(wt)-vs-epo-with-hmgb1-(eh))
        - [Untreated (WT) vs HMGB1 (H) ](#untreated-(wt)-vs-hmgb1-(h))
        - [EPO with Vehicle (EV) vs EPO with HMGB1 (EH) ](#epo-with-vehicle-(ev)-vs-epo-with-hmgb1-(eh))
        - [EPO with Vehicle (EV) vs HMGB1 (H) ](#epo-with-vehicle-(ev)-vs-hmgb1-(h))
        - [EPO with HMGB1 (EH) vs HMGB1 (H) ](#epo-with-hmgb1-(eh)-vs-hmgb1-(h))


```r
# Load libraries
set.seed(0)
source("https://bioconductor.org/biocLite.R")
library("DESeq2")
library("RColorBrewer")
library("vsn")
library("pheatmap")
library("EnhancedVolcano")
library("gridExtra")
library("grid")
library("ggplot2")
library("gplots")
library("affy")
library("reshape2")
library("gProfileR")
library("ComplexHeatmap")
```
***
## Summary of Experiment:

1. Flow-sorted Colony-Forming Units of Erythroid (CFU-E) cells. The cells are then starved for six hours to stop the Erythropoietin (EPO). Then these cells are treated with:
    - EPO + Vehicle (Peripheral Blood Smear (PBS))
    - EPO + High mobility group box 1 (HMGB1) protein
    - HMGB1

2. The treatment occurs for sixty minutes and then RNA was extracted. Three biological replicates are used for reproducibility.

3. Single-end stranded RNA-seq performed on wild-type EPO starved cells and treated with 1,2,3 as stated above.

4. Labels:
    - WT => Untreated
    - EV => EPO + Vehicle
    - EH => EPO + HMGB1
    - H => HMGB1
    - E1 => Samples from Patient 1
    - E2 => Samples from Patient 2
    - E3 => Samples from Patient 3

## Exploring data with DESeq2

<!-- ### Basic DESeq2 setup -->


```r
directory <- "/Users/rshah22/Documents/Projects/RNAseq_Erythrocid/counts_analysis/"
files <- read.table(paste(directory, 
    "deseq2_sampleinfo.txt", sep = ""), 
    header = T)
condition <- read.csv(paste(directory, 
    "deseq2_phenotype.txt", sep = ""), 
    header = T)
sampleTable <- data.frame(sampleName = files$sampleName, 
    fileName = files$fileName, 
    condition.type = condition$condition, 
    condition.treatment = condition$treatment, 
    condition.experiment = condition$experiment)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
    directory = directory, design = ~condition.type + 
        condition.experiment)
# ddsHTSeq
ddsHTSeq$condition.type <- factor(ddsHTSeq$condition.type, 
    levels = c("WT", "EV", "EH", 
        "H"))
ddsHTSeq$condition.experiment <- factor(ddsHTSeq$condition.experiment, 
    levels = c("E1", "E2", "E3"))
```


### Running DESeq2


```r
dds <- DESeq(ddsHTSeq)
```


### Make counts table to explore that data

```r
count.table <- counts(dds)
count.table <- count.table[-c((grep("^HIST|RNA|^RP11", 
    rownames(count.table)))), ]
stats.per.sample <- data.frame(t(do.call(cbind, 
    lapply(count.table, summary))))
stats.per.sample$libsum <- apply(count.table, 
    2, sum)  ## libsum
stats.per.sample$perc05 <- apply(count.table, 
    2, quantile, 0.05)
stats.per.sample$perc10 <- apply(count.table, 
    2, quantile, 0.1)
stats.per.sample$perc90 <- apply(count.table, 
    2, quantile, 0.9)
stats.per.sample$perc95 <- apply(count.table, 
    2, quantile, 0.95)
stats.per.sample$zeros <- apply(count.table == 
    0, 2, sum)
stats.per.sample$percent.zeros <- 100 * 
    stats.per.sample$zeros/nrow(count.table)
kable(stats.per.sample[sample(1:ncol(count.table), 
    size = 10), ], caption = "**Table: statistics per sample. ** We only display a random selection of 10 samples. ")
```



Table: **Table: statistics per sample. ** We only display a random selection of 10 samples. 

      Min.   X1st.Qu.   Median   Mean   X3rd.Qu.   Max.     libsum   perc05   perc10   perc90   perc95   zeros   percent.zeros
---  -----  ---------  -------  -----  ---------  -----  ---------  -------  -------  -------  -------  ------  --------------
11       0          0        0      0          0      0   16025018        0        0    864.0   1727.2   24354        54.92929
3        5          5        5      5          5      5   15796424        0        0    800.4   1621.0   24715        55.74351
4        0          0        0      0          0      0   15273891        0        0    774.0   1569.2   24872        56.09762
6       11         11       11     11         11     11   13073035        0        0    645.0   1292.0   24544        55.35783
8        0          0        0      0          0      0   16006088        0        0    804.0   1620.0   24365        54.95410
2        3          3        3      3          3      3   15509458        0        0    793.0   1603.0   24712        55.73674
9        0          0        0      0          0      0   14905408        0        0    789.0   1569.0   24223        54.63383
5        1          1        1      1          1      1   16355877        0        0    811.4   1628.2   24136        54.43760
12     882        882      882    882        882    882   13745586        0        0    728.0   1477.0   24754        55.83147
7       33         33       33     33         33     33   16044909        0        0    801.0   1620.2   24460        55.16837

***
<!--### Assign color based on sample type -->

```r
col.sampletype <- c(WT = "#9E0142", 
    EV = "#D0384D", EH = "#EE6445", 
    H = "#FA9C58")
expDesign <- condition
expDesign$color <- col.sampletype[as.vector(expDesign$condition)]
```

### Histograms of counts per gene {.tabset .tabset-fade .tabset-pills}
#### matrix

```r
plot.new()
hist(as.matrix(count.table), col = "blue", 
    border = "white", breaks = 100)
```

![](deseq2_analysis_files/figure-html/code block 6-1.png)<!-- -->

```r
invisible(dev.off())
```

#### truncated matrix

```r
plot.new()
hist(as.matrix(count.table), col = "blue", 
    border = "white", breaks = 20000, 
    xlim = c(0, 500), main = "Counts per gene", 
    xlab = "Counts (truncated axis)", 
    ylab = "Number of genes", las = 1, 
    cex.axis = 0.7)
```

![](deseq2_analysis_files/figure-html/code block 7-1.png)<!-- -->

```r
invisible(dev.off())
```

#### log2-transformation

```r
epsilon <- 1  # pseudo-count to avoid problems with log(0)
plot.new()
hist(as.matrix(log2(count.table + 
    epsilon)), breaks = 100, col = "blue", 
    border = "white", main = "Log2-transformed counts per gene", 
    xlab = "log2(counts+1)", ylab = "Number of genes", 
    las = 1, cex.axis = 0.7)
```

![](deseq2_analysis_files/figure-html/code block 8-1.png)<!-- -->

```r
invisible(dev.off())
```

#### box-plot

```r
plot.new()
boxplot(log2(count.table + epsilon), 
    col = expDesign$color, pch = ".", 
    horizontal = TRUE, cex.axis = 0.5, 
    las = 1, ylab = "Samples", 
    xlab = "log2(Counts +1)")
```

![](deseq2_analysis_files/figure-html/code block 9-1.png)<!-- -->

```r
invisible(dev.off())
```
#### density-plots

```r
plot.new()
plotDensity(log2(count.table + 
    epsilon), lty = 1, col = expDesign$color, 
    lwd = 2)
grid()
legend("topright", legend = names(col.sampletype), 
    col = col.sampletype, lwd = 2)
```

![](deseq2_analysis_files/figure-html/code block 10-1.png)<!-- -->

```r
invisible(dev.off())
```



#### scatter-plots
-  all points are aligned along the diagonal, with a relatively wider dispersion at the bottom, corresponding to small number fluctuations.
- comparing across samples types there are points discarding from the diagonal.

```r
plotFun <- function(x, y) {
    dns <- densCols(x, y)
    points(x, y, col = dns, pch = ".", 
        panel.first = grid())
    abline(a = 0, b = 1, col = "brown")
}
plot.new()
pairs(log2(count.table[, sample(ncol(count.table), 
    12)] + epsilon), panel = plotFun, 
    lower.panel = NULL)
```

![](deseq2_analysis_files/figure-html/code block 11-1.png)<!-- -->

```r
invisible(dev.off())
```

### Eliminate undetected genes

```r
prop.null <- apply(count.table, 
    2, function(x) 100 * mean(x == 
        0))
plot.new()
barplot(prop.null, main = "Percentage of null counts per sample", 
    horiz = TRUE, cex.names = 0.5, 
    las = 1, col = expDesign$color, 
    ylab = "Samples", xlab = "% of null counts")
```

![](deseq2_analysis_files/figure-html/code block 12-1.png)<!-- -->

```r
invisible(dev.off())
count.table <- count.table[rowSums(count.table) > 
    0, ]
```

###  Re-run DESeq2 & Normalization 

```r
dds0 <- DESeqDataSetFromMatrix(countData = count.table, 
    colData = expDesign, design = ~condition + 
        experiment)
dds.norm <- estimateSizeFactors(dds0)
# sizeFactors(dds.norm)
```

#### Plot normalized data {.tabset .tabset-fade .tabset-pills}
##### box-plot

```r
par(mfrow = c(1, 2), cex.lab = 0.7)
boxplot(log2(counts(dds.norm) + 
    epsilon), col = col.sampletype, 
    cex.axis = 0.7, las = 1, xlab = "log2(counts)", 
    horizontal = TRUE, main = "Raw counts")
boxplot(log2(counts(dds.norm, normalized = TRUE) + 
    epsilon), col = col.sampletype, 
    cex.axis = 0.7, las = 1, xlab = "log2(normalized counts)", 
    horizontal = TRUE, main = "Normalized counts")
```

![](deseq2_analysis_files/figure-html/code block 14-1.png)<!-- -->

```r
invisible(dev.off())
```

##### density-plot

```r
par(mfrow = c(1, 2), cex.lab = 0.7)
plotDensity(log2(counts(dds.norm) + 
    epsilon), col = col.sampletype, 
    xlab = "log2(counts)", cex.lab = 0.7)
plotDensity(log2(counts(dds.norm, 
    normalized = TRUE) + epsilon), 
    col = col.sampletype, xlab = "log2(normalized counts)", 
    cex.lab = 0.7)
```

![](deseq2_analysis_files/figure-html/code block 15-1.png)<!-- -->

```r
invisible(dev.off())
```

### Count variance is related to mean

```r
## Computing mean and variance
norm.counts <- counts(dds.norm, 
    normalized = TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 
    1, var)
## sum(mean.counts==0) # Number
## of completely undetected
## genes
norm.counts.stats <- data.frame(min = apply(norm.counts, 
    2, min), mean = apply(norm.counts, 
    2, mean), median = apply(norm.counts, 
    2, median), max = apply(norm.counts, 
    2, max), zeros = apply(norm.counts == 
    0, 2, sum), percent.zeros = 100 * 
    apply(norm.counts == 0, 2, 
        sum)/nrow(norm.counts), 
    perc05 = apply(norm.counts, 
        2, quantile, 0.05), perc10 = apply(norm.counts, 
        2, quantile, 0.1), perc90 = apply(norm.counts, 
        2, quantile, 0.9), perc95 = apply(norm.counts, 
        2, quantile, 0.95))
kable(norm.counts.stats)
```

       min       mean     median        max   zeros   percent.zeros   perc05   perc10     perc90     perc95
----  ----  ---------  ---------  ---------  ------  --------------  -------  -------  ---------  ---------
S01      0   573.3060   10.42074   302131.5    6882        26.11269        0        0   1358.864   2394.581
S02      0   568.1465   10.61987   302188.5    6730        25.53595        0        0   1333.856   2346.992
S03      0   579.9790   10.64411   321006.9    6733        25.54733        0        0   1346.963   2362.604
S04      0   581.0388   10.02579   350783.2    6890        26.14305        0        0   1359.497   2363.379
S05      0   572.3780   11.98990   374807.9    6154        23.35041        0        0   1294.909   2275.314
S06      0   568.7111   12.61162   383578.0    6562        24.89850        0        0   1278.360   2241.773
S07      0   579.0621   12.36500   380738.5    6478        24.57978        0        0   1304.032   2315.395
S08      0   584.3166   11.54536   366004.1    6383        24.21931        0        0   1328.870   2354.868
S09      0   543.7968   12.49969   260979.1    6241        23.68052        0        0   1306.698   2255.713
S10      0   560.4016   11.20687   265525.5    6653        25.24379        0        0   1335.656   2337.041
S11      0   555.9397   10.97169   230312.2    6372        24.17758        0        0   1360.489   2350.684
S12      0   567.6020   10.88288   264424.5    6772        25.69531        0        0   1381.255   2399.239

```r
## Mean and variance
## relationship
mean.var.col <- densCols(x = log2(mean.counts), 
    y = log2(variance.counts))
plot.new()
plot(x = log2(mean.counts), y = log2(variance.counts), 
    pch = 16, cex = 0.5, col = mean.var.col, 
    main = "Mean-variance relationship", 
    xlab = "Mean log2(normalized counts) per gene", 
    ylab = "Variance of log2(normalized counts)", 
    panel.first = grid())
abline(a = 0, b = 1, col = "brown")
```

![](deseq2_analysis_files/figure-html/code block 16-1.png)<!-- -->

```r
invisible(dev.off())
```

### Estimated Dispersion for each gene

- Shows the mean of normalized counts (x axis) and dispersion estimate for each genes

```r
## Performing estimation of
## dispersion parameter
dds.disp <- estimateDispersions(dds.norm)
plotDispEsts(dds.disp)
```

![](deseq2_analysis_files/figure-html/code block 17-1.png)<!-- -->

```r
invisible(dev.off())
```

### PCA and heatmap on all samples using reguralized log transformation {.tabset .tabset-fade .tabset-pills}
 
 - Transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.

#### PCA using rlog

```r
rld <- rlog(dds.disp, blind = TRUE)
hmcol <- colorRampPalette(brewer.pal(11, 
    "Spectral"))(12)
plotPCA(rld, ntop = 5000, intgroup = c("condition", 
    "experiment")) + scale_color_manual(values = hmcol)
```

![](deseq2_analysis_files/figure-html/code block 18-1.png)<!-- -->

```r
invisible(dev.off())
```

#### Heatmap of count matrix for first 50 genes

```r
select <- order(rowMedians(counts(dds.disp, 
    normalized = TRUE)), decreasing = TRUE)[1:50]
df <- as.data.frame(colData(dds.disp)[, 
    c("condition", "experiment")])
# Heatmap of count matrix
pheatmap(assay(rld)[select, ], 
    cluster_rows = TRUE, show_rownames = TRUE, 
    cluster_cols = TRUE, annotation_col = df)
```

![](deseq2_analysis_files/figure-html/code block 19-1.png)<!-- -->

```r
invisible(dev.off())
```

#### Heatmap of sample distance based on rlog

```r
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, 
    rld$experiment, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, 
    "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, 
    clustering_distance_cols = sampleDists, 
    col = colors)
```

![](deseq2_analysis_files/figure-html/code block 20-1.png)<!-- -->

```r
invisible(dev.off())
```

### Perform differential expression call

#### Set standard cutoff and Enhanced Volcano Plots

```r
p_cutoff = 0.1
fc_cutoff = 0.5
xlim <- c(1, 1e+05)
ylim <- c(-6, 6)
wald.test <- nbinomWaldTest(dds.disp)
enhanced_volcano_plots <- function(title, 
    results, pvalue, foldchange) {
    p = EnhancedVolcano(results, 
        lab = rownames(results), 
        x = "log2FoldChange", y = "padj", 
        xlab = bquote(~Log[2] ~ 
            "fold change"), ylab = bquote(~-Log[10] ~ 
            adjusted ~ italic(P)), 
        pCutoff = pvalue, FCcutoff = foldchange, 
        xlim = c(-6, 6), transcriptLabSize = 3, 
        title = title, colAlpha = 1, 
        legend = c("NS", "Log2 FC", 
            "Adjusted p-value", 
            "Adjusted p-value & Log2 FC"), 
        legendPosition = "bottom", 
        legendLabSize = 10, legendIconSize = 3, 
        DrawConnectors = TRUE, 
        widthConnectors = 0.5, 
        colConnectors = "black")
    return(p)
}
plotCounts_gg <- function(i, dds, 
    intgroup) {
    group <- if (length(intgroup) == 
        1) {
        colData(dds)[[intgroup]]
    } else if (length(intgroup) == 
        2) {
        lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]), 
            levels(colData(dds)[[intgroup[2]]]), 
            function(x, y) paste(x, 
                y, sep = " : "))))
        droplevels(factor(apply(as.data.frame(colData(dds)[, 
            intgroup, drop = FALSE]), 
            1, paste, collapse = " : "), 
            levels = lvls))
    } else {
        factor(apply(as.data.frame(colData(dds)[, 
            intgroup, drop = FALSE]), 
            1, paste, collapse = " : "))
    }
    data <- plotCounts(dds, gene = i, 
        intgroup = intgroup, returnData = TRUE)
    data <- cbind(data, data.frame(group = group))
    main <- rownames(dds)[i]
    ggplot(data, aes(x = group, 
        y = count)) + geom_boxplot() + 
        ylab("Normalized count") + 
        ggtitle(main) + coord_trans(y = "log2") + 
        scale_x_discrete(limits = c("WT", 
            "EV", "EH", "H"))
}
```

#### Untreated (WT) vs EPO with Vehicle (EV) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_WT_EV <- results(wald.test, 
    contrast = c("condition", "WT", 
        "EV"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_WT_EV_ashr <- lfcShrink(wald.test, 
    contrast = c("condition", "WT", 
        "EV"), res = res_WT_EV, 
    type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_WT_EV, xlim = xlim, 
    ylim = ylim, main = "normal")
plotMA(res_WT_EV_ashr, xlim = xlim, 
    ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 22-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_WT_EV$padj, breaks = 20, 
    col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 23-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_WT_EV <- res_WT_EV[order(res_WT_EV$padj), 
    ]
write.csv(as.data.frame(resOrdered_WT_EV), 
    file = "condition_WT_EV_alpha0.1_results.csv")
WT_EV_p1 <- enhanced_volcano_plots("WT vs. EV (padj=0.05,log2fc=0.5)", 
    res_WT_EV, 0.05, fc_cutoff)
WT_EV_p2 <- enhanced_volcano_plots("WT vs. EV (padj=0.1,log2fc=0.5)", 
    res_WT_EV, p_cutoff, fc_cutoff)
grid.arrange(WT_EV_p1, WT_EV_p2, 
    nrow = 2, ncol = 1)
grid.rect(gp = gpar(fill = NA))
```

![](deseq2_analysis_files/figure-html/code block 24-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Plot counts

```r
nBestFeatures = 20
ord <- order(res_WT_EV$padj, decreasing = FALSE)
for (i in head(ord, nBestFeatures)) {
    print(plotCounts_gg(i, dds = wald.test, 
        intgroup = c("condition")))
}
```

![](deseq2_analysis_files/figure-html/code block 25-1.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-2.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-3.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-4.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-5.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-6.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-7.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-8.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-9.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-10.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-11.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-12.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-13.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-14.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-15.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-16.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-17.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-18.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-19.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 25-20.png)<!-- -->

##### Select genes based on FDR and make heatmap


```r
  gene.kept <- rownames(res_WT_EV)[res_WT_EV$padj <= p_cutoff & !is.na(res_WT_EV$padj) & ( res_WT_EV$log2FoldChange <= -0.5 | res_WT_EV$log2FoldChange >= 0.5)]
  count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
  heatmap.2(as.matrix(count.table.kept), 
            scale="row", 
            hclust=function(x) hclust(x,method="average"), 
            distfun=function(x) as.dist((1-cor(t(x)))/2), 
            trace="none", 
            density="none", 
            #labRow="",
            cexCol=0.7)
```

![](deseq2_analysis_files/figure-html/code black 26-1.png)<!-- -->

```r
  invisible(dev.off())
```

##### Do functional enrichment


```r
res_WT_EV.df <- na.omit(data.frame(res_WT_EV))
induced.sign <- rownames(res_WT_EV.df)[res_WT_EV.df$log2FoldChange >= 
    0.5 & res_WT_EV.df$padj < p_cutoff]
# head(induced.sign)
# names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, 
        organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), 
        ]
    # term.induced$p.value
    kable(term.induced[1:10, c("term.name", 
        "term.size", "query.size", 
        "overlap.size", "recall", 
        "precision", "p.value", 
        "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis wit gProfileR. ** ")
}
```

```
## No genes found that have induced expression
```

```r
repressed.sign <- rownames(res_WT_EV.df)[res_WT_EV.df$log2FoldChange <= 
    -0.5 & res_WT_EV.df$padj < 
    p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, 
        organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), 
        ]
    kable(term.repressed[1:10, 
        c("term.name", "term.size", 
            "query.size", "overlap.size", 
            "recall", "precision", 
            "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, 
            digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ")
}
```



Table: **Table: repressed genes functional analysis with gProfileR. ** 

     term.name                                                            term.size   query.size   overlap.size   recall   precision    p.value  intersection                                                          
---  ------------------------------------------------------------------  ----------  -----------  -------------  -------  ----------  ---------  ----------------------------------------------------------------------
12   Jak-STAT signaling pathway                                                 162           22              8    0.049       0.364   7.00e-07  IL4R,OSM,CISH,SOCS2,CDKN1A,MYC,PIM1,SOCS3                             
17   Signaling by Interleukins                                                  463           21             10    0.022       0.476   1.60e-06  IL4R,OSM,CISH,SOCS2,CD80,CDKN1A,MYC,PIM1,SOCS3,IRS2                   
16   Cytokine Signaling in Immune system                                        685           21             11    0.016       0.524   4.60e-06  IL4R,OSM,CISH,EGR1,SOCS2,CD80,CDKN1A,MYC,PIM1,SOCS3,IRS2              
18   Interleukin-4 and 13 signaling                                             111           21              6    0.054       0.286   1.68e-05  IL4R,OSM,CDKN1A,MYC,PIM1,SOCS3                                        
20   Growth hormone receptor signaling                                           24           21              4    0.167       0.190   3.57e-05  CISH,SOCS2,SOCS3,IRS2                                                 
8    response to insulin                                                        255           37              7    0.027       0.189   1.06e-03  CISH,EGR1,SOCS2,MYC,SOCS3,IRS2,PPARA                                  
15   Immune System                                                             2010           21             13    0.006       0.619   5.50e-03  IL4R,OSM,CISH,EGR1,SOCS2,CD80,CDKN1A,MYC,PIM1,MB21D1,HSPA6,SOCS3,IRS2 
10   response to oxygen-containing compound                                    1555           37             13    0.008       0.351   7.17e-03  CPEB4,CISH,KLF9,EGR1,SOCS2,CDKN1A,MYC,PIM1,AQP3,RARG,SOCS3,IRS2,PPARA 
14   TFAP2 (AP-2) family regulates transcription of cell cycle factors            5           21              2    0.400       0.095   1.11e-02  CDKN1A,MYC                                                            
19   Interleukin-7 signaling                                                     36           21              3    0.083       0.143   1.38e-02  CISH,SOCS2,IRS2                                                       

```r
# kable(head(term.induced[,c('p.value',
# 'term.name','intersection')],
# 10))
```

#### Untreated (WT) vs EPO with HMGB1 (EH)

##### Plots {.tabset .tabset-fade .tabset-pills}
###### MA-plot

```r
res_WT_EH <- results(wald.test, 
    contrast = c("condition", "WT", 
        "EH"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_WT_EH_ashr <- lfcShrink(wald.test, 
    contrast = c("condition", "WT", 
        "EH"), res = res_WT_EV, 
    type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_WT_EH, xlim = xlim, 
    ylim = ylim, main = "normal")
plotMA(res_WT_EH_ashr, xlim = xlim, 
    ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 28-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_WT_EH$padj, breaks = 20, 
    col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 29-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_WT_EH <- res_WT_EH[order(res_WT_EH$padj), 
    ]
write.csv(as.data.frame(resOrdered_WT_EH), 
    file = "condition_WT_EH_alpha0.1_results.csv")
WT_EH_p1 <- enhanced_volcano_plots("WT vs. EH (padj=0.05,log2fc=0.5)", 
    res_WT_EH, 0.05, fc_cutoff)
WT_EH_p2 <- enhanced_volcano_plots("WT vs. EH (padj=0.1,log2fc=0.5)", 
    res_WT_EH, p_cutoff, fc_cutoff)
grid.arrange(WT_EH_p1, WT_EH_p2, 
    nrow = 2, ncol = 1)
grid.rect(gp = gpar(fill = NA))
```

![](deseq2_analysis_files/figure-html/code block 30-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Plot counts

```r
nBestFeatures = 20
ord <- order(res_WT_EH$padj, decreasing = FALSE)
for (i in head(ord, nBestFeatures)) {
    print(plotCounts_gg(i, dds = wald.test, 
        intgroup = c("condition")))
}
```

![](deseq2_analysis_files/figure-html/code block 31-1.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-2.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-3.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-4.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-5.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-6.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-7.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-8.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-9.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-10.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-11.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-12.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-13.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-14.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-15.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-16.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-17.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-18.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-19.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 31-20.png)<!-- -->

##### Select genes based on FDR and make heatmap


```r
  gene.kept <- rownames(res_WT_EH)[res_WT_EH$padj <= p_cutoff & !is.na(res_WT_EH$padj) & ( res_WT_EH$log2FoldChange <= -0.5 | res_WT_EH$log2FoldChange >= 0.5)]
  count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
  heatmap.2(as.matrix(count.table.kept), 
            scale="row", 
            hclust=function(x) hclust(x,method="average"), 
            distfun=function(x) as.dist((1-cor(t(x)))/2), 
            trace="none", 
            density="none", 
            #labRow="",
            cexCol=0.7)
```

![](deseq2_analysis_files/figure-html/code black 32-1.png)<!-- -->

```r
  invisible(dev.off())
```

##### Do functional enrichment


```r
res_WT_EH.df <- na.omit(data.frame(res_WT_EH))
induced.sign <- rownames(res_WT_EV.df)[res_WT_EH.df$log2FoldChange >= 
    0.5 & res_WT_EH.df$padj < p_cutoff]
# head(induced.sign)
# names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, 
        organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), 
        ]
    # term.induced$p.value
    kable(term.induced[1:10, c("term.name", 
        "term.size", "query.size", 
        "overlap.size", "recall", 
        "precision", "p.value", 
        "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis wit gProfileR. ** ")
}
```



Table: **Table: induced gene functional analysis wit gProfileR. ** 

       term.name                                        term.size   query.size   overlap.size   recall   precision   p.value  intersection                                                                                                                                                                                                                     
-----  ----------------------------------------------  ----------  -----------  -------------  -------  ----------  --------  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1      macrolide binding                                       16           45              3    0.188       0.067    0.0201  FKBP4,FKBP3,FKBP2                                                                                                                                                                                                                
2      FK506 binding                                           16           45              3    0.188       0.067    0.0201  FKBP4,FKBP3,FKBP2                                                                                                                                                                                                                
3      Factor: E2F-3; motif: GGCGGGN; match class: 1         9297           45             35    0.004       0.778    0.0464  FKBP4,RANBP3,SNX24,XPO1,ZFAND6,KIF16B,FKBP3,TRIB3,RELB,FAM184A,ZNF142,LGALSL,ACTR10,TET1,COL6A1,DCAF6,CDCA7,TAGLN,HNRNPDL,TRIM11,OTUD6B,EYA3,ERCC3,DACT1,FAM69B,SPATA33,DDB1,UPF3A,DUS1L,ZNF524,FKBP2,TALDO1,LIN28B,DACT3,ZNF521 
NA     NA                                                      NA           NA             NA       NA          NA        NA  NA                                                                                                                                                                                                                               
NA.1   NA                                                      NA           NA             NA       NA          NA        NA  NA                                                                                                                                                                                                                               
NA.2   NA                                                      NA           NA             NA       NA          NA        NA  NA                                                                                                                                                                                                                               
NA.3   NA                                                      NA           NA             NA       NA          NA        NA  NA                                                                                                                                                                                                                               
NA.4   NA                                                      NA           NA             NA       NA          NA        NA  NA                                                                                                                                                                                                                               
NA.5   NA                                                      NA           NA             NA       NA          NA        NA  NA                                                                                                                                                                                                                               
NA.6   NA                                                      NA           NA             NA       NA          NA        NA  NA                                                                                                                                                                                                                               

```r
repressed.sign <- rownames(res_WT_EH.df)[res_WT_EH.df$log2FoldChange <= 
    -0.5 & res_WT_EH.df$padj < 
    p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, 
        organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), 
        ]
    kable(term.repressed[1:10, 
        c("term.name", "term.size", 
            "query.size", "overlap.size", 
            "recall", "precision", 
            "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, 
            digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ")
}
```



Table: **Table: repressed genes functional analysis with gProfileR. ** 

     term.name                                                                                                                    term.size   query.size   overlap.size   recall   precision    p.value  intersection                                                                                                                                                                                                                        
---  --------------------------------------------------------------------------------------------------------------------------  ----------  -----------  -------------  -------  ----------  ---------  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
19   Cytokine Signaling in Immune system                                                                                                685           22             11    0.016       0.500   7.90e-06  ICAM1,OSM,NFKBIA,RELB,CISH,EGR1,SOCS2,IRF1,MYC,PIM1,JUN                                                                                                                                                                             
18   Immune System                                                                                                                     2010           22             16    0.008       0.727   1.80e-05  ICAM1,OSM,NFKBIA,RELB,CISH,TNFAIP3,EGR1,SOCS2,IRF1,MYC,PIM1,NFKBIE,MB21D1,SLCO4C1,JUN,UBOX5                                                                                                                                         
16   Epstein-Barr virus infection                                                                                                       197           26              7    0.036       0.269   2.44e-04  ICAM1,NFKBIA,RELB,TNFAIP3,MYC,NFKBIE,JUN                                                                                                                                                                                            
20   Signaling by Interleukins                                                                                                          463           22              8    0.017       0.364   6.51e-04  ICAM1,OSM,NFKBIA,CISH,SOCS2,MYC,PIM1,JUN                                                                                                                                                                                            
3    I-kappaB/NF-kappaB complex                                                                                                           7           42              3    0.429       0.071   8.70e-04  BCL3,NFKBIA,RELB                                                                                                                                                                                                                    
22   Factor: E2F; motif: GGCGSG; match class: 1                                                                                       10205           42             37    0.004       0.881   1.14e-03  MSMO1,VMP1,KLF6,BCL3,ICAM1,OSM,NFKBIA,RELB,ARRDC3,CISH,TNFAIP3,EGR1,SOCS2,ARL4A,IRF1,A4GALT,TMEM160,NINJ1,BHLHE40,STX11,MYC,PIM1,NFKBIZ,NFKBIE,RAB11FIP1,ETS2,ANKRD33B,MB21D1,KCNK5,LETM2,STARD5,SLCO4C1,JUN,DDN,H1FX,UBOX5,TSC22D2 
15   HTLV-I infection                                                                                                                   253           26              7    0.028       0.269   1.26e-03  ICAM1,NFKBIA,RELB,EGR1,MYC,ETS2,JUN                                                                                                                                                                                                 
14   C-type lectin receptor signaling pathway                                                                                           104           26              5    0.048       0.192   1.70e-03  BCL3,NFKBIA,RELB,IRF1,JUN                                                                                                                                                                                                           
13   TNF signaling pathway                                                                                                              108           26              5    0.046       0.192   2.04e-03  BCL3,ICAM1,NFKBIA,TNFAIP3,JUN                                                                                                                                                                                                       
9    TNF-alpha/NF-kappa B signaling complex (CHUK, KPNA3, NFKB2, NFKBIB, REL, IKBKG, NFKB1, NFKBIE, RELB, NFKBIA, RELA, TNIP2)           12           12              3    0.250       0.250   2.54e-03  NFKBIA,RELB,NFKBIE                                                                                                                                                                                                                  

```r
# kable(head(term.induced[,c('p.value',
# 'term.name','intersection')],
# 10))
```

***

#### Untreated (WT) vs HMGB1 (H) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_WT_H <- results(wald.test, 
    contrast = c("condition", "WT", 
        "H"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_WT_H_ashr <- lfcShrink(wald.test, 
    contrast = c("condition", "WT", 
        "H"), res = res_WT_H, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_WT_H, xlim = xlim, ylim = ylim, 
    main = "normal")
plotMA(res_WT_H_ashr, xlim = xlim, 
    ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 34-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_WT_H$padj, breaks = 20, 
    col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 35-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_WT_H <- res_WT_H[order(res_WT_H$padj), 
    ]
write.csv(as.data.frame(resOrdered_WT_H), 
    file = "condition_WT_H_alpha0.1_results.csv")
WT_H_p1 <- enhanced_volcano_plots("WT vs. H (padj=0.05,log2fc=0.5)", 
    res_WT_H, 0.05, fc_cutoff)
WT_H_p2 <- enhanced_volcano_plots("WT vs. H (padj=0.1,log2fc=0.5)", 
    res_WT_H, p_cutoff, fc_cutoff)
grid.arrange(WT_H_p1, WT_H_p2, 
    nrow = 2, ncol = 1)
grid.rect(gp = gpar(fill = NA))
```

![](deseq2_analysis_files/figure-html/code block 36-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Plot counts

```r
nBestFeatures = 20
ord <- order(res_WT_H$padj, decreasing = FALSE)
for (i in head(ord, nBestFeatures)) {
    print(plotCounts_gg(i, dds = wald.test, 
        intgroup = c("condition")))
}
```

![](deseq2_analysis_files/figure-html/code block 37-1.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-2.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-3.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-4.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-5.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-6.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-7.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-8.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-9.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-10.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-11.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-12.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-13.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-14.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-15.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-16.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-17.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-18.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-19.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 37-20.png)<!-- -->

##### Select genes based on FDR and make heatmap


```r
  gene.kept <- rownames(res_WT_H)[res_WT_H$padj <= p_cutoff & !is.na(res_WT_H$padj) & ( res_WT_H$log2FoldChange <= -0.5 | res_WT_H$log2FoldChange >= 0.5)]
  count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
  heatmap.2(as.matrix(count.table.kept), 
            scale="row", 
            hclust=function(x) hclust(x,method="average"), 
            distfun=function(x) as.dist((1-cor(t(x)))/2), 
            trace="none", 
            density="none", 
            #labRow="",
            cexCol=0.7)
```

![](deseq2_analysis_files/figure-html/code black 38-1.png)<!-- -->

```r
  invisible(dev.off())
```

##### Do functional enrichment


```r
res_WT_H.df <- na.omit(data.frame(res_WT_H))
induced.sign <- rownames(res_WT_H.df)[res_WT_H.df$log2FoldChange >= 
    0.5 & res_WT_H.df$padj < p_cutoff]
# head(induced.sign)
# names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, 
        organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), 
        ]
    # term.induced$p.value
    kable(term.induced[1:10, c("term.name", 
        "term.size", "query.size", 
        "overlap.size", "recall", 
        "precision", "p.value", 
        "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ")
}
```



Table: **Table: induced gene functional analysis with gProfileR. ** 

       term.name   term.size   query.size   overlap.size   recall   precision   p.value   intersection 
-----  ----------  ----------  -----------  -------------  -------  ----------  --------  -------------
NA     NA          NA          NA           NA             NA       NA          NA        NA           
NA.1   NA          NA          NA           NA             NA       NA          NA        NA           
NA.2   NA          NA          NA           NA             NA       NA          NA        NA           
NA.3   NA          NA          NA           NA             NA       NA          NA        NA           
NA.4   NA          NA          NA           NA             NA       NA          NA        NA           
NA.5   NA          NA          NA           NA             NA       NA          NA        NA           
NA.6   NA          NA          NA           NA             NA       NA          NA        NA           
NA.7   NA          NA          NA           NA             NA       NA          NA        NA           
NA.8   NA          NA          NA           NA             NA       NA          NA        NA           
NA.9   NA          NA          NA           NA             NA       NA          NA        NA           

```r
repressed.sign <- rownames(res_WT_H.df)[res_WT_H.df$log2FoldChange <= 
    -0.5 & res_WT_H.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, 
        organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), 
        ]
    kable(term.repressed[1:10, 
        c("term.name", "term.size", 
            "query.size", "overlap.size", 
            "recall", "precision", 
            "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, 
            digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ")
}
```



Table: **Table: repressed genes functional analysis with gProfileR. ** 

     term.name                                                                                                                    term.size   query.size   overlap.size   recall   precision    p.value  intersection                                                                                                    
---  --------------------------------------------------------------------------------------------------------------------------  ----------  -----------  -------------  -------  ----------  ---------  ----------------------------------------------------------------------------------------------------------------
11   TNF signaling pathway                                                                                                              108           14              5    0.046       0.357   6.81e-05  BCL3,ICAM1,NFKBIA,TNFAIP3,JUN                                                                                   
7    Epstein-Barr virus infection                                                                                                       197           14              5    0.025       0.357   1.31e-03  ICAM1,NFKBIA,TNFAIP3,NFKBIE,JUN                                                                                 
3    TNF-alpha/NF-kappa B signaling complex (CHUK, KPNA3, NFKB2, NFKBIB, REL, IKBKG, NFKB1, NFKBIE, RELB, NFKBIA, RELA, TNIP2)           12            6              2    0.167       0.333   1.11e-02  NFKBIA,NFKBIE                                                                                                   
4    CHUK-NFKB2-REL-IKBKG-SPAG9-NFKB1-NFKBIE-COPB2-TNIP1-NFKBIA-RELA-TNIP2 complex                                                       12            6              2    0.167       0.333   1.11e-02  NFKBIA,NFKBIE                                                                                                   
2    regulation of DNA binding transcription factor activity                                                                            400           23              6    0.015       0.261   1.34e-02  ICAM1,NFKBIA,TNFAIP3,BHLHE40,NFKBIE,JUN                                                                         
6    B cell receptor signaling pathway                                                                                                   70           14              3    0.043       0.214   1.68e-02  NFKBIA,NFKBIE,JUN                                                                                               
1    nucleotide-binding oligomerization domain containing 1 signaling pathway                                                             4           23              2    0.500       0.087   1.93e-02  NFKBIA,TNFAIP3                                                                                                  
12   Factor: NF-kappaB; motif: GGGGATYCCC                                                                                              6390           22             17    0.003       0.773   3.27e-02  BCL3,ICAM1,HIVEP1,NFKBIA,ARRDC3,TNFAIP3,NINJ1,STX11,NFKBIZ,NFKBIE,ZFP36L2,KLF15,TPSAB1,SLCO4C1,FAM174A,JUN,H1F0 
13   Factor: NF-kappaB; motif: GGGGATYCCC; match class: 0                                                                              6390           22             17    0.003       0.773   3.27e-02  BCL3,ICAM1,HIVEP1,NFKBIA,ARRDC3,TNFAIP3,NINJ1,STX11,NFKBIZ,NFKBIE,ZFP36L2,KLF15,TPSAB1,SLCO4C1,FAM174A,JUN,H1F0 
8    Th1 and Th2 cell differentiation                                                                                                    90           14              3    0.033       0.214   3.53e-02  NFKBIA,NFKBIE,JUN                                                                                               

```r
# kable(head(term.induced[,c('p.value',
# 'term.name','intersection')],
# 10))
```

#### EPO with Vehicle (EV) vs EPO with HMGB1 (EH) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_EV_EH <- results(wald.test, 
    contrast = c("condition", "EV", 
        "EH"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_EV_EH_ashr <- lfcShrink(wald.test, 
    contrast = c("condition", "EV", 
        "EH"), res = res_EV_EH, 
    type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_EV_EH, xlim = xlim, 
    ylim = ylim, main = "normal")
plotMA(res_EV_EH_ashr, xlim = xlim, 
    ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 40-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_EV_EH$padj, breaks = 20, 
    col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 41-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_EV_EH <- res_EV_EH[order(res_EV_EH$padj), 
    ]
write.csv(as.data.frame(resOrdered_EV_EH), 
    file = "condition_EV_EH_alpha0.1_results.csv")
EV_EH_p1 <- enhanced_volcano_plots("EV vs. EH (padj=0.05,log2fc=0.5)", 
    res_EV_EH, 0.05, fc_cutoff)
EV_EH_p2 <- enhanced_volcano_plots("EV vs. EH (padj=0.1,log2fc=0.5)", 
    res_EV_EH, p_cutoff, fc_cutoff)
grid.arrange(EV_EH_p1, EV_EH_p2, 
    nrow = 2, ncol = 1)
grid.rect(gp = gpar(fill = NA))
```

![](deseq2_analysis_files/figure-html/code block 42-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Plot counts

```r
nBestFeatures = 20
ord <- order(res_EV_EH$padj, decreasing = FALSE)
for (i in head(ord, nBestFeatures)) {
    print(plotCounts_gg(i, dds = wald.test, 
        intgroup = c("condition")))
}
```

![](deseq2_analysis_files/figure-html/code block 43-1.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-2.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-3.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-4.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-5.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-6.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-7.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-8.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-9.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-10.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-11.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-12.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-13.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-14.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-15.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-16.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-17.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-18.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-19.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 43-20.png)<!-- -->

##### Select genes based on FDR and make heatmap


```r
  gene.kept <- rownames(res_EV_EH)[res_EV_EH$padj <= p_cutoff & !is.na(res_EV_EH$padj) & ( res_EV_EH$log2FoldChange <= -0.5 | res_EV_EH$log2FoldChange >= 0.5)]
  count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
  heatmap.2(as.matrix(count.table.kept), 
            scale="row", 
            hclust=function(x) hclust(x,method="average"), 
            distfun=function(x) as.dist((1-cor(t(x)))/2), 
            trace="none", 
            density="none", 
            #labRow="",
            cexCol=0.7)
```

![](deseq2_analysis_files/figure-html/code black 44-1.png)<!-- -->

```r
  invisible(dev.off())
```

##### Do functional enrichment


```r
res_EV_EH.df <- na.omit(data.frame(res_EV_EH))
induced.sign <- rownames(res_EV_EH.df)[res_EV_EH.df$log2FoldChange >= 
    0.5 & res_EV_EH.df$padj < p_cutoff]
# head(induced.sign)
# names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, 
        organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), 
        ]
    # term.induced$p.value
    kable(term.induced[1:10, c("term.name", 
        "term.size", "query.size", 
        "overlap.size", "recall", 
        "precision", "p.value", 
        "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ")
}
```



Table: **Table: induced gene functional analysis with gProfileR. ** 

       term.name                                             term.size   query.size   overlap.size   recall   precision   p.value  intersection                                                             
-----  ---------------------------------------------------  ----------  -----------  -------------  -------  ----------  --------  -------------------------------------------------------------------------
3      positive regulation of tau-protein kinase activity            4           25              2    0.500       0.080    0.0200  EGR1,NAB2                                                                
1      regulation of cellular protein metabolic process           2438           25             12    0.005       0.480    0.0344  GADD45B,EIF5,EEF2K,GPLD1,EGR1,ID1,HSPA2,ARRDC4,DUSP2,NAB2,SH3BP5L,HSPA1A 
4      Attenuation phase                                            14           15              2    0.143       0.133    0.0439  HSPA2,HSPA1A                                                             
2      regulation of tau-protein kinase activity                     6           25              2    0.333       0.080    0.0499  EGR1,NAB2                                                                
NA     NA                                                           NA           NA             NA       NA          NA        NA  NA                                                                       
NA.1   NA                                                           NA           NA             NA       NA          NA        NA  NA                                                                       
NA.2   NA                                                           NA           NA             NA       NA          NA        NA  NA                                                                       
NA.3   NA                                                           NA           NA             NA       NA          NA        NA  NA                                                                       
NA.4   NA                                                           NA           NA             NA       NA          NA        NA  NA                                                                       
NA.5   NA                                                           NA           NA             NA       NA          NA        NA  NA                                                                       

```r
repressed.sign <- rownames(res_EV_EH.df)[res_EV_EH.df$log2FoldChange <= 
    -0.5 & res_EV_EH.df$padj < 
    p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, 
        organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), 
        ]
    kable(term.repressed[1:10, 
        c("term.name", "term.size", 
            "query.size", "overlap.size", 
            "recall", "precision", 
            "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, 
            digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ")
}
```



Table: **Table: repressed genes functional analysis with gProfileR. ** 

     term.name                                                                                                                    term.size   query.size   overlap.size   recall   precision    p.value  intersection                                  
---  --------------------------------------------------------------------------------------------------------------------------  ----------  -----------  -------------  -------  ----------  ---------  ----------------------------------------------
28   TNF signaling pathway                                                                                                              108           15              7    0.065       0.467   0.00e+00  BIRC3,TRAF1,BCL3,ICAM1,NFKBIA,TNFAIP3,JUN     
31   NF-kappa B signaling pathway                                                                                                        93           15              6    0.065       0.400   9.00e-07  BIRC3,TRAF1,ICAM1,NFKBIA,RELB,TNFAIP3         
30   Epstein-Barr virus infection                                                                                                       197           15              7    0.036       0.467   2.70e-06  TRAF1,ICAM1,NFKBIA,RELB,TNFAIP3,NFKBIE,JUN    
21   I-kappaB/NF-kappaB complex                                                                                                           7           20              3    0.429       0.150   7.49e-05  BCL3,NFKBIA,RELB                              
1    regulation of DNA binding transcription factor activity                                                                            400           20              7    0.018       0.350   2.14e-04  TRAF1,ICAM1,NFKBIA,TNFAIP3,BHLHE40,NFKBIE,JUN 
23   TNF-alpha/NF-kappa B signaling complex (CHUK, KPNA3, NFKB2, NFKBIB, REL, IKBKG, NFKB1, NFKBIE, RELB, NFKBIA, RELA, TNIP2)           12            8              3    0.250       0.375   3.12e-04  NFKBIA,RELB,NFKBIE                            
41   TNFR1-induced NFkappaB signaling pathway                                                                                            26           12              3    0.115       0.250   7.68e-04  BIRC3,TRAF1,TNFAIP3                           
15   cellular response to tumor necrosis factor                                                                                         290           20              6    0.021       0.300   7.94e-04  BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3      
14   response to tumor necrosis factor                                                                                                  307           20              6    0.020       0.300   1.11e-03  BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3      
40   Regulation of TNFR1 signaling                                                                                                       33           12              3    0.091       0.250   1.61e-03  BIRC3,TRAF1,TNFAIP3                           

```r
# kable(head(term.induced[,c('p.value',
# 'term.name','intersection')],
# 10))
```

#### EPO with Vehicle (EV) vs HMGB1 (H) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_EV_H <- results(wald.test, 
    contrast = c("condition", "EV", 
        "H"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_EV_H_ashr <- lfcShrink(wald.test, 
    contrast = c("condition", "EV", 
        "H"), res = res_EV_H, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_EV_H, xlim = xlim, ylim = ylim, 
    main = "normal")
plotMA(res_EV_H_ashr, xlim = xlim, 
    ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 46-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_EV_H$padj, breaks = 20, 
    col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 47-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_EV_H <- res_EV_H[order(res_EV_H$padj), 
    ]
write.csv(as.data.frame(resOrdered_EV_H), 
    file = "condition_EV_H_alpha0.1_results.csv")
EV_H_p1 <- enhanced_volcano_plots("EV vs. H (padj=0.05,log2fc=0.5)", 
    res_EV_H, 0.05, fc_cutoff)
EV_H_p2 <- enhanced_volcano_plots("EV vs. H (padj=0.1,log2fc=0.5)", 
    res_EV_H, p_cutoff, fc_cutoff)
grid.arrange(EV_H_p1, EV_H_p2, 
    nrow = 2, ncol = 1)
grid.rect(gp = gpar(fill = NA))
```

![](deseq2_analysis_files/figure-html/code block 48-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Plot counts

```r
nBestFeatures = 20
ord <- order(res_EV_H$padj, decreasing = FALSE)
for (i in head(ord, nBestFeatures)) {
    print(plotCounts_gg(i, dds = wald.test, 
        intgroup = c("condition")))
}
```

![](deseq2_analysis_files/figure-html/code block 49-1.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-2.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-3.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-4.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-5.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-6.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-7.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-8.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-9.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-10.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-11.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-12.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-13.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-14.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-15.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-16.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-17.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-18.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-19.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 49-20.png)<!-- -->

##### Select genes based on FDR and make heatmap


```r
  gene.kept <- rownames(res_EV_H)[res_EV_H$padj <= p_cutoff & !is.na(res_EV_H$padj) & ( res_EV_H$log2FoldChange <= -0.5 | res_EV_H$log2FoldChange >= 0.5)]
  count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
  heatmap.2(as.matrix(count.table.kept), 
            scale="row", 
            hclust=function(x) hclust(x,method="average"), 
            distfun=function(x) as.dist((1-cor(t(x)))/2), 
            trace="none", 
            density="none", 
            #labRow="",
            cexCol=0.7)
```

![](deseq2_analysis_files/figure-html/code black 50-1.png)<!-- -->

```r
  invisible(dev.off())
```

##### Do functional enrichment


```r
res_EV_H.df <- na.omit(data.frame(res_EV_H))
induced.sign <- rownames(res_EV_H.df)[res_EV_H.df$log2FoldChange >= 
    0.5 & res_EV_H.df$padj < p_cutoff]
# head(induced.sign)
# names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, 
        organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), 
        ]
    # term.induced$p.value
    kable(term.induced[1:10, c("term.name", 
        "term.size", "query.size", 
        "overlap.size", "recall", 
        "precision", "p.value", 
        "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ")
}
```



Table: **Table: induced gene functional analysis with gProfileR. ** 

     term.name                                      term.size   query.size   overlap.size   recall   precision    p.value  intersection                                                     
---  --------------------------------------------  ----------  -----------  -------------  -------  ----------  ---------  -----------------------------------------------------------------
36   Growth hormone receptor signaling                     24           35              5    0.208       0.143   3.80e-06  CISH,SOCS2,SOCS3,SOCS1,IRS2                                      
35   Cytokine Signaling in Immune system                  685           35             12    0.018       0.343   3.36e-04  OSM,CISH,EGR1,SOCS2,CD80,CDKN1A,PIM1,EIF4A1,FOS,SOCS3,SOCS1,IRS2 
38   Interleukin-4 and 13 signaling                       111           35              6    0.054       0.171   4.71e-04  OSM,CDKN1A,PIM1,FOS,SOCS3,SOCS1                                  
37   Signaling by Interleukins                            463           35             10    0.022       0.286   5.17e-04  OSM,CISH,SOCS2,CD80,CDKN1A,PIM1,FOS,SOCS3,SOCS1,IRS2             
33   Jak-STAT signaling pathway                           162           35              7    0.043       0.200   6.45e-04  OSM,CISH,SOCS2,CDKN1A,PIM1,SOCS3,SOCS1                           
25   negative regulation of kinase activity               272           66              9    0.033       0.136   1.01e-03  GADD45B,CISH,SOCS2,CDKN1A,DUSP2,SH3BP5L,SOCS3,SOCS1,IRS2         
32   Prolactin signaling pathway                           70           35              5    0.071       0.143   1.27e-03  CISH,SOCS2,FOS,SOCS3,SOCS1                                       
39   Interleukin-7 signaling                               36           35              4    0.111       0.114   1.74e-03  CISH,SOCS2,SOCS1,IRS2                                            
24   negative regulation of transferase activity          298           66              9    0.030       0.136   2.18e-03  GADD45B,CISH,SOCS2,CDKN1A,DUSP2,SH3BP5L,SOCS3,SOCS1,IRS2         
31   protein kinase inhibitor activity                     95           66              6    0.063       0.091   2.52e-03  CISH,SOCS2,CDKN1A,SH3BP5L,SOCS3,SOCS1                            

```r
repressed.sign <- rownames(res_EV_H.df)[res_EV_H.df$log2FoldChange <= 
    -0.5 & res_EV_H.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, 
        organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), 
        ]
    kable(term.repressed[1:10, 
        c("term.name", "term.size", 
            "query.size", "overlap.size", 
            "recall", "precision", 
            "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, 
            digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ")
}
```



Table: **Table: repressed genes functional analysis with gProfileR. ** 

     term.name                                           term.size   query.size   overlap.size   recall   precision    p.value  intersection                                     
---  -------------------------------------------------  ----------  -----------  -------------  -------  ----------  ---------  -------------------------------------------------
17   TNF signaling pathway                                     108           17              7    0.065       0.412   1.00e-07  BIRC3,TRAF1,BCL3,ICAM1,NFKBIA,TNFAIP3,JUN        
15   NF-kappa B signaling pathway                               93           17              5    0.054       0.294   9.21e-05  BIRC3,TRAF1,ICAM1,NFKBIA,TNFAIP3                 
3    cellular response to tumor necrosis factor                290           24              7    0.024       0.292   9.55e-05  BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3,ZFP36L2 
2    response to tumor necrosis factor                         307           24              7    0.023       0.292   1.41e-04  BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3,ZFP36L2 
10   liver development                                         139           24              5    0.036       0.208   1.47e-03  NFKBIA,TNFAIP3,SLCO2B1,CITED2,JUN                
21   TNFR1-induced NFkappaB signaling pathway                   26           15              3    0.115       0.200   1.53e-03  BIRC3,TRAF1,TNFAIP3                              
9    hepaticobiliary system development                        142           24              5    0.035       0.208   1.64e-03  NFKBIA,TNFAIP3,SLCO2B1,CITED2,JUN                
20   Regulation of TNFR1 signaling                              33           15              3    0.091       0.200   3.20e-03  BIRC3,TRAF1,TNFAIP3                              
16   Epstein-Barr virus infection                              197           17              5    0.025       0.294   3.60e-03  TRAF1,ICAM1,NFKBIA,TNFAIP3,JUN                   
4    tumor necrosis factor-mediated signaling pathway          174           24              5    0.029       0.208   4.46e-03  BIRC3,TNFRSF9,TRAF1,NFKBIA,TNFAIP3               

```r
# kable(head(term.induced[,c('p.value',
# 'term.name','intersection')],
# 10))
```

#### EPO with HMGB1 (EH) vs HMGB1 (H) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_EH_H <- results(wald.test, 
    contrast = c("condition", "EH", 
        "H"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_EH_H_ashr <- lfcShrink(wald.test, 
    contrast = c("condition", "EH", 
        "H"), res = res_EH_H, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_EH_H, xlim = xlim, ylim = ylim, 
    main = "normal")
plotMA(res_EH_H_ashr, xlim = xlim, 
    ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 52-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_EH_H$padj, breaks = 20, 
    col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 53-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_EH_H <- res_EH_H[order(res_EH_H$padj), 
    ]
write.csv(as.data.frame(resOrdered_EH_H), 
    file = "condition_EH_H_alpha0.1_results.csv")
EH_H_p1 <- enhanced_volcano_plots("EH vs. H (padj=0.05,log2fc=0.5)", 
    res_EH_H, 0.05, fc_cutoff)
EH_H_p2 <- enhanced_volcano_plots("EH vs. H (padj=0.1,log2fc=0.5)", 
    res_EH_H, p_cutoff, fc_cutoff)
grid.arrange(EH_H_p1, EH_H_p2, 
    nrow = 2, ncol = 1)
grid.rect(gp = gpar(fill = NA))
```

![](deseq2_analysis_files/figure-html/code block 54-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Plot counts

```r
nBestFeatures = 20
ord <- order(res_EH_H$padj, decreasing = FALSE)
for (i in head(ord, nBestFeatures)) {
    print(plotCounts_gg(i, dds = wald.test, 
        intgroup = c("condition")))
}
```

![](deseq2_analysis_files/figure-html/code block 55-1.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-2.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-3.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-4.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-5.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-6.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-7.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-8.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-9.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-10.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-11.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-12.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-13.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-14.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-15.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-16.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-17.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-18.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-19.png)<!-- -->![](deseq2_analysis_files/figure-html/code block 55-20.png)<!-- -->

##### Select genes based on FDR and make heatmap


```r
  gene.kept <- rownames(res_EH_H)[res_EH_H$padj <= p_cutoff & !is.na(res_EH_H$padj) & ( res_EH_H$log2FoldChange <= -0.5 | res_EH_H$log2FoldChange >= 0.5)]
  count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
  heatmap.2(as.matrix(count.table.kept), 
            scale="row", 
            hclust=function(x) hclust(x,method="average"), 
            distfun=function(x) as.dist((1-cor(t(x)))/2), 
            trace="none", 
            density="none", 
            #labRow="",
            cexCol=0.7)
```

![](deseq2_analysis_files/figure-html/code black 56-1.png)<!-- -->

```r
  invisible(dev.off())
```

##### Do functional enrichment


```r
res_EH_H.df <- na.omit(data.frame(res_EH_H))
induced.sign <- rownames(res_EH_H.df)[res_EH_H.df$log2FoldChange >= 
    0.5 & res_EH_H.df$padj < p_cutoff]
# head(induced.sign)
# names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, 
        organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), 
        ]
    # term.induced$p.value
    kable(term.induced[1:10, c("term.name", 
        "term.size", "query.size", 
        "overlap.size", "recall", 
        "precision", "p.value", 
        "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ")
}
```



Table: **Table: induced gene functional analysis with gProfileR. ** 

       term.name                              term.size   query.size   overlap.size   recall   precision    p.value  intersection             
-----  ------------------------------------  ----------  -----------  -------------  -------  ----------  ---------  -------------------------
3      Jak-STAT signaling pathway                   162            6              4    0.025       0.667   0.000146  OSM,CISH,SOCS2,PIM1      
4      Cytokine Signaling in Immune system          685            7              5    0.007       0.714   0.004300  OSM,CISH,EGR1,SOCS2,PIM1 
5      Growth hormone receptor signaling             24            7              2    0.083       0.286   0.020700  CISH,SOCS2               
6      Signaling by Interleukins                    463            7              4    0.009       0.571   0.023100  OSM,CISH,SOCS2,PIM1      
7      Interleukin-7 signaling                       36            7              2    0.056       0.286   0.047100  CISH,SOCS2               
1      EGR-EP300 complex                              2            2              1    0.500       0.500   0.049800  EGR1                     
2      RCP-Rab11 complex                              2            2              1    0.500       0.500   0.049800  RAB11FIP1                
NA     NA                                            NA           NA             NA       NA          NA         NA  NA                       
NA.1   NA                                            NA           NA             NA       NA          NA         NA  NA                       
NA.2   NA                                            NA           NA             NA       NA          NA         NA  NA                       

```r
repressed.sign <- rownames(res_EH_H.df)[res_EH_H.df$log2FoldChange <= 
    -0.5 & res_EH_H.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, 
        organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), 
        ]
    kable(term.repressed[1:10, 
        c("term.name", "term.size", 
            "query.size", "overlap.size", 
            "recall", "precision", 
            "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, 
            digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ")
}
```

```
## No genes found that have repressed expression
```

```r
# kable(head(term.induced[,c('p.value',
# 'term.name','intersection')],
# 10))
```

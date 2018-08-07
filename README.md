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
library("kableExtra")
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
files <- read.table(paste(directory, "deseq2_sampleinfo.txt", sep = ""), header = T)
condition <- read.csv(paste(directory, "deseq2_phenotype.txt", sep = ""), header = T)
sampleTable <- data.frame(sampleName = files$sampleName, fileName = files$fileName, 
    condition.type = condition$condition, condition.treatment = condition$treatment, 
    condition.experiment = condition$experiment)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, 
    design = ~condition.type + condition.experiment)
# ddsHTSeq
ddsHTSeq$condition.type <- factor(ddsHTSeq$condition.type, levels = c("WT", 
    "EV", "EH", "H"))
ddsHTSeq$condition.experiment <- factor(ddsHTSeq$condition.experiment, levels = c("E1", 
    "E2", "E3"))
```


### Running DESeq2


```r
dds <- DESeq(ddsHTSeq)
```


### Make counts table to explore that data

```r
count.table <- counts(dds)
count.table <- count.table[-c((grep("^HIST|RNA|^RP11|^LINC|^RPL|^RNVU", rownames(count.table)))), 
    ]
stats.per.sample <- data.frame(t(do.call(cbind, lapply(count.table, summary))))
stats.per.sample$libsum <- apply(count.table, 2, sum)  ## libsum
stats.per.sample$perc05 <- apply(count.table, 2, quantile, 0.05)
stats.per.sample$perc10 <- apply(count.table, 2, quantile, 0.1)
stats.per.sample$perc90 <- apply(count.table, 2, quantile, 0.9)
stats.per.sample$perc95 <- apply(count.table, 2, quantile, 0.95)
stats.per.sample$zeros <- apply(count.table == 0, 2, sum)
stats.per.sample$percent.zeros <- 100 * stats.per.sample$zeros/nrow(count.table)
kable(stats.per.sample[sample(1:ncol(count.table), size = 10), ], caption = "**Table: statistics per sample. ** We only display a random selection of 10 samples. ") %>% 
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: statistics per sample. ** We only display a random selection of 10 samples. </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Min. </th>
   <th style="text-align:right;"> X1st.Qu. </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> X3rd.Qu. </th>
   <th style="text-align:right;"> Max. </th>
   <th style="text-align:right;"> libsum </th>
   <th style="text-align:right;"> perc05 </th>
   <th style="text-align:right;"> perc10 </th>
   <th style="text-align:right;"> perc90 </th>
   <th style="text-align:right;"> perc95 </th>
   <th style="text-align:right;"> zeros </th>
   <th style="text-align:right;"> percent.zeros </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 11 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 15453771 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 884 </td>
   <td style="text-align:right;"> 1739.5 </td>
   <td style="text-align:right;"> 23465 </td>
   <td style="text-align:right;"> 54.47981 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 15188212 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 818 </td>
   <td style="text-align:right;"> 1632.0 </td>
   <td style="text-align:right;"> 23824 </td>
   <td style="text-align:right;"> 55.31332 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 14681229 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 793 </td>
   <td style="text-align:right;"> 1579.0 </td>
   <td style="text-align:right;"> 23971 </td>
   <td style="text-align:right;"> 55.65462 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 12569247 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 661 </td>
   <td style="text-align:right;"> 1301.0 </td>
   <td style="text-align:right;"> 23649 </td>
   <td style="text-align:right;"> 54.90701 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 15417872 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 824 </td>
   <td style="text-align:right;"> 1626.5 </td>
   <td style="text-align:right;"> 23469 </td>
   <td style="text-align:right;"> 54.48910 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 14920843 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 812 </td>
   <td style="text-align:right;"> 1617.0 </td>
   <td style="text-align:right;"> 23813 </td>
   <td style="text-align:right;"> 55.28778 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 14310464 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 808 </td>
   <td style="text-align:right;"> 1580.0 </td>
   <td style="text-align:right;"> 23343 </td>
   <td style="text-align:right;"> 54.19656 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 15712110 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 828 </td>
   <td style="text-align:right;"> 1648.0 </td>
   <td style="text-align:right;"> 23228 </td>
   <td style="text-align:right;"> 53.92956 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:right;"> 882 </td>
   <td style="text-align:right;"> 882 </td>
   <td style="text-align:right;"> 882 </td>
   <td style="text-align:right;"> 882 </td>
   <td style="text-align:right;"> 882 </td>
   <td style="text-align:right;"> 882 </td>
   <td style="text-align:right;"> 13250203 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 749 </td>
   <td style="text-align:right;"> 1487.5 </td>
   <td style="text-align:right;"> 23863 </td>
   <td style="text-align:right;"> 55.40387 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 15438596 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 819 </td>
   <td style="text-align:right;"> 1631.0 </td>
   <td style="text-align:right;"> 23545 </td>
   <td style="text-align:right;"> 54.66555 </td>
  </tr>
</tbody>
</table>

***
<!--### Assign color based on sample type -->

```r
col.sampletype <- c(WT = "#9E0142", EV = "#D0384D", EH = "#EE6445", H = "#FA9C58")
expDesign <- condition
expDesign$color <- col.sampletype[as.vector(expDesign$condition)]
```

### Histograms of counts per gene {.tabset .tabset-fade .tabset-pills}
#### matrix

```r
plot.new()
hist(as.matrix(count.table), col = "blue", border = "white", breaks = 100)
```

![](deseq2_analysis_files/figure-html/code block 6-1.png)<!-- -->

```r
invisible(dev.off())
```

#### truncated matrix

```r
plot.new()
hist(as.matrix(count.table), col = "blue", border = "white", breaks = 20000, 
    xlim = c(0, 500), main = "Counts per gene", xlab = "Counts (truncated axis)", 
    ylab = "Number of genes", las = 1, cex.axis = 0.7)
```

![](deseq2_analysis_files/figure-html/code block 7-1.png)<!-- -->

```r
invisible(dev.off())
```

#### log2-transformation

```r
epsilon <- 1  # pseudo-count to avoid problems with log(0)
plot.new()
hist(as.matrix(log2(count.table + epsilon)), breaks = 100, col = "blue", border = "white", 
    main = "Log2-transformed counts per gene", xlab = "log2(counts+1)", ylab = "Number of genes", 
    las = 1, cex.axis = 0.7)
```

![](deseq2_analysis_files/figure-html/code block 8-1.png)<!-- -->

```r
invisible(dev.off())
```

#### box-plot

```r
plot.new()
boxplot(log2(count.table + epsilon), col = expDesign$color, pch = ".", horizontal = TRUE, 
    cex.axis = 0.5, las = 1, ylab = "Samples", xlab = "log2(Counts +1)")
```

![](deseq2_analysis_files/figure-html/code block 9-1.png)<!-- -->

```r
invisible(dev.off())
```
#### density-plots

```r
plot.new()
plotDensity(log2(count.table + epsilon), lty = 1, col = expDesign$color, lwd = 2)
grid()
legend("topright", legend = names(col.sampletype), col = col.sampletype, lwd = 2)
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
    points(x, y, col = dns, pch = ".", panel.first = grid())
    abline(a = 0, b = 1, col = "brown")
}
plot.new()
pairs(log2(count.table[, sample(ncol(count.table), 12)] + epsilon), panel = plotFun, 
    lower.panel = NULL)
```

![](deseq2_analysis_files/figure-html/code block 11-1.png)<!-- -->

```r
invisible(dev.off())
```

### Eliminate undetected genes

```r
prop.null <- apply(count.table, 2, function(x) 100 * mean(x == 0))
plot.new()
barplot(prop.null, main = "Percentage of null counts per sample", horiz = TRUE, 
    cex.names = 0.5, las = 1, col = expDesign$color, ylab = "Samples", xlab = "% of null counts")
```

![](deseq2_analysis_files/figure-html/code block 12-1.png)<!-- -->

```r
invisible(dev.off())
count.table <- count.table[rowSums(count.table) > 0, ]
```

###  Re-run DESeq2 & Normalization 

```r
dds0 <- DESeqDataSetFromMatrix(countData = count.table, colData = expDesign, 
    design = ~condition + experiment)
dds.norm <- estimateSizeFactors(dds0)
# sizeFactors(dds.norm)
```

#### Plot normalized data {.tabset .tabset-fade .tabset-pills}
##### box-plot

```r
par(mfrow = c(1, 2), cex.lab = 0.7)
boxplot(log2(counts(dds.norm) + epsilon), col = col.sampletype, cex.axis = 0.7, 
    las = 1, xlab = "log2(counts)", horizontal = TRUE, main = "Raw counts")
boxplot(log2(counts(dds.norm, normalized = TRUE) + epsilon), col = col.sampletype, 
    cex.axis = 0.7, las = 1, xlab = "log2(normalized counts)", horizontal = TRUE, 
    main = "Normalized counts")
```

![](deseq2_analysis_files/figure-html/code block 14-1.png)<!-- -->

```r
invisible(dev.off())
```

##### density-plot

```r
par(mfrow = c(1, 2), cex.lab = 0.7)
plotDensity(log2(counts(dds.norm) + epsilon), col = col.sampletype, xlab = "log2(counts)", 
    cex.lab = 0.7)
plotDensity(log2(counts(dds.norm, normalized = TRUE) + epsilon), col = col.sampletype, 
    xlab = "log2(normalized counts)", cex.lab = 0.7)
```

![](deseq2_analysis_files/figure-html/code block 15-1.png)<!-- -->

```r
invisible(dev.off())
```

### Count variance is related to mean

```r
## Computing mean and variance
norm.counts <- counts(dds.norm, normalized = TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)
## sum(mean.counts==0) # Number of completely undetected genes
norm.counts.stats <- data.frame(min = apply(norm.counts, 2, min), mean = apply(norm.counts, 
    2, mean), median = apply(norm.counts, 2, median), max = apply(norm.counts, 
    2, max), zeros = apply(norm.counts == 0, 2, sum), percent.zeros = 100 * 
    apply(norm.counts == 0, 2, sum)/nrow(norm.counts), perc05 = apply(norm.counts, 
    2, quantile, 0.05), perc10 = apply(norm.counts, 2, quantile, 0.1), perc90 = apply(norm.counts, 
    2, quantile, 0.9), perc95 = apply(norm.counts, 2, quantile, 0.95))
kable(norm.counts.stats, caption = "**Table: statistics for count. ** ") %>% 
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: statistics for count. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> min </th>
   <th style="text-align:right;"> mean </th>
   <th style="text-align:right;"> median </th>
   <th style="text-align:right;"> max </th>
   <th style="text-align:right;"> zeros </th>
   <th style="text-align:right;"> percent.zeros </th>
   <th style="text-align:right;"> perc05 </th>
   <th style="text-align:right;"> perc10 </th>
   <th style="text-align:right;"> perc90 </th>
   <th style="text-align:right;"> perc95 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> S01 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 563.8231 </td>
   <td style="text-align:right;"> 11.46740 </td>
   <td style="text-align:right;"> 302252.5 </td>
   <td style="text-align:right;"> 6643 </td>
   <td style="text-align:right;"> 25.80106 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1365.663 </td>
   <td style="text-align:right;"> 2374.794 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S02 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 559.4646 </td>
   <td style="text-align:right;"> 11.58476 </td>
   <td style="text-align:right;"> 302174.0 </td>
   <td style="text-align:right;"> 6489 </td>
   <td style="text-align:right;"> 25.20294 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1339.392 </td>
   <td style="text-align:right;"> 2330.178 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S03 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 570.9703 </td>
   <td style="text-align:right;"> 11.61488 </td>
   <td style="text-align:right;"> 321093.4 </td>
   <td style="text-align:right;"> 6500 </td>
   <td style="text-align:right;"> 25.24566 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1355.457 </td>
   <td style="text-align:right;"> 2351.433 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S04 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 571.9695 </td>
   <td style="text-align:right;"> 11.03392 </td>
   <td style="text-align:right;"> 350959.9 </td>
   <td style="text-align:right;"> 6647 </td>
   <td style="text-align:right;"> 25.81660 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1367.203 </td>
   <td style="text-align:right;"> 2347.918 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S05 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 563.0195 </td>
   <td style="text-align:right;"> 12.91646 </td>
   <td style="text-align:right;"> 374931.7 </td>
   <td style="text-align:right;"> 5904 </td>
   <td style="text-align:right;"> 22.93083 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1302.164 </td>
   <td style="text-align:right;"> 2257.336 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S06 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 560.1683 </td>
   <td style="text-align:right;"> 13.76947 </td>
   <td style="text-align:right;"> 383893.9 </td>
   <td style="text-align:right;"> 6325 </td>
   <td style="text-align:right;"> 24.56597 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1289.740 </td>
   <td style="text-align:right;"> 2226.867 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S07 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 570.1387 </td>
   <td style="text-align:right;"> 13.31151 </td>
   <td style="text-align:right;"> 380605.6 </td>
   <td style="text-align:right;"> 6221 </td>
   <td style="text-align:right;"> 24.16204 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1311.945 </td>
   <td style="text-align:right;"> 2298.233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S08 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 576.0588 </td>
   <td style="text-align:right;"> 12.50583 </td>
   <td style="text-align:right;"> 365956.1 </td>
   <td style="text-align:right;"> 6145 </td>
   <td style="text-align:right;"> 23.86686 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1336.584 </td>
   <td style="text-align:right;"> 2334.453 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S09 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 534.5886 </td>
   <td style="text-align:right;"> 13.46544 </td>
   <td style="text-align:right;"> 261061.3 </td>
   <td style="text-align:right;"> 6019 </td>
   <td style="text-align:right;"> 23.37748 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1314.227 </td>
   <td style="text-align:right;"> 2242.669 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 551.5089 </td>
   <td style="text-align:right;"> 12.22801 </td>
   <td style="text-align:right;"> 265576.0 </td>
   <td style="text-align:right;"> 6403 </td>
   <td style="text-align:right;"> 24.86892 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1341.005 </td>
   <td style="text-align:right;"> 2321.997 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S11 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 548.5404 </td>
   <td style="text-align:right;"> 11.88076 </td>
   <td style="text-align:right;"> 230210.7 </td>
   <td style="text-align:right;"> 6141 </td>
   <td style="text-align:right;"> 23.85132 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1367.567 </td>
   <td style="text-align:right;"> 2335.665 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S12 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 559.3352 </td>
   <td style="text-align:right;"> 11.95553 </td>
   <td style="text-align:right;"> 264079.3 </td>
   <td style="text-align:right;"> 6539 </td>
   <td style="text-align:right;"> 25.39713 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1390.102 </td>
   <td style="text-align:right;"> 2375.673 </td>
  </tr>
</tbody>
</table>

```r
## Mean and variance relationship
mean.var.col <- densCols(x = log2(mean.counts), y = log2(variance.counts))
plot.new()
plot(x = log2(mean.counts), y = log2(variance.counts), pch = 16, cex = 0.5, 
    col = mean.var.col, main = "Mean-variance relationship", xlab = "Mean log2(normalized counts) per gene", 
    ylab = "Variance of log2(normalized counts)", panel.first = grid())
abline(a = 0, b = 1, col = "brown")
```

![](deseq2_analysis_files/figure-html/code block 16-1.png)<!-- -->

```r
invisible(dev.off())
```

### Estimated Dispersion for each gene

- Shows the mean of normalized counts (x axis) and dispersion estimate for each genes

```r
## Performing estimation of dispersion parameter
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
hmcol <- colorRampPalette(brewer.pal(11, "Spectral"))(12)
plotPCA(rld, ntop = 5000, intgroup = c("condition", "experiment")) + scale_color_manual(values = hmcol)
```

![](deseq2_analysis_files/figure-html/code block 18-1.png)<!-- -->

```r
invisible(dev.off())
```

#### Heatmap of count matrix for first 50 genes

```r
select <- order(rowMedians(counts(dds.disp, normalized = TRUE)), decreasing = TRUE)[1:50]
df <- as.data.frame(colData(dds.disp)[, c("condition", "experiment")])
# Heatmap of count matrix
pheatmap(assay(rld)[select, ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
    annotation_col = df)
```

![](deseq2_analysis_files/figure-html/code block 19-1.png)<!-- -->

```r
invisible(dev.off())
```

#### Heatmap of sample distance based on rlog

```r
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$experiment, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, 
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
enhanced_volcano_plots <- function(title, results, pvalue, foldchange) {
    p = EnhancedVolcano(results, lab = rownames(results), x = "log2FoldChange", 
        y = "padj", xlab = bquote(~Log[2] ~ "fold change"), ylab = bquote(~-Log[10] ~ 
            adjusted ~ italic(P)), pCutoff = pvalue, FCcutoff = foldchange, 
        xlim = c(-6, 6), transcriptLabSize = 3, title = title, colAlpha = 1, 
        legend = c("NS", "Log2 FC", "Adjusted p-value", "Adjusted p-value & Log2 FC"), 
        legendPosition = "bottom", legendLabSize = 10, legendIconSize = 3, DrawConnectors = TRUE, 
        widthConnectors = 0.5, colConnectors = "black")
    return(p)
}
plotCounts_gg <- function(i, dds, intgroup) {
    group <- if (length(intgroup) == 1) {
        colData(dds)[[intgroup]]
    } else if (length(intgroup) == 2) {
        lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]), levels(colData(dds)[[intgroup[2]]]), 
            function(x, y) paste(x, y, sep = " : "))))
        droplevels(factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]), 
            1, paste, collapse = " : "), levels = lvls))
    } else {
        factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]), 
            1, paste, collapse = " : "))
    }
    data <- plotCounts(dds, gene = i, intgroup = intgroup, returnData = TRUE)
    data <- cbind(data, data.frame(group = group))
    main <- rownames(dds)[i]
    ggplot(data, aes(x = group, y = count)) + geom_boxplot() + ylab("Normalized count") + 
        ggtitle(main) + coord_trans(y = "log2") + scale_x_discrete(limits = c("WT", 
        "EV", "EH", "H"))
}
# for Kable
rowlim = 10
```

#### Untreated (WT) vs EPO with Vehicle (EV) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_WT_EV <- results(wald.test, contrast = c("condition", "WT", "EV"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_WT_EV_ashr <- lfcShrink(wald.test, contrast = c("condition", "WT", "EV"), 
    res = res_WT_EV, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_WT_EV, xlim = xlim, ylim = ylim, main = "normal")
plotMA(res_WT_EV_ashr, xlim = xlim, ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 22-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_WT_EV$padj, breaks = 20, col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 23-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_WT_EV <- res_WT_EV[order(res_WT_EV$padj), ]
write.csv(as.data.frame(resOrdered_WT_EV), file = "condition_WT_EV_alpha0.1_results.csv")
WT_EV_p1 <- enhanced_volcano_plots("WT vs. EV (padj=0.05,log2fc=0.5)", res_WT_EV, 
    0.05, fc_cutoff)
WT_EV_p2 <- enhanced_volcano_plots("WT vs. EV (padj=0.1,log2fc=0.5)", res_WT_EV, 
    p_cutoff, fc_cutoff)
grid.arrange(WT_EV_p1, WT_EV_p2, nrow = 2, ncol = 1)
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
    print(plotCounts_gg(i, dds = wald.test, intgroup = c("condition")))
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
induced.sign <- rownames(res_WT_EV.df)[res_WT_EV.df$log2FoldChange >= 0.5 & 
    res_WT_EV.df$padj < p_cutoff]
# head(induced.sign) names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), ]
    na.omit(term.induced)
    # term.induced$p.value
    if (nrow(term.induced) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.induced)
    }
    kable(term.induced[1:rowlim, c("term.name", "term.size", "query.size", "overlap.size", 
        "recall", "precision", "p.value", "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis wit gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

```
## No genes found that have induced expression
```

```r
repressed.sign <- rownames(res_WT_EV.df)[res_WT_EV.df$log2FoldChange <= -0.5 & 
    res_WT_EV.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), ]
    na.omit(term.repressed)
    if (nrow(term.repressed) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.repressed)
    }
    kable(term.repressed[1:rowlim, c("term.name", "term.size", "query.size", 
        "overlap.size", "recall", "precision", "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: repressed genes functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> Jak-STAT signaling pathway </td>
   <td style="text-align:right;"> 162 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.049 </td>
   <td style="text-align:right;"> 0.364 </td>
   <td style="text-align:right;"> 7.00e-07 </td>
   <td style="text-align:left;"> IL4R,OSM,CISH,SOCS2,CDKN1A,MYC,PIM1,SOCS3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:left;"> Signaling by Interleukins </td>
   <td style="text-align:right;"> 463 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.022 </td>
   <td style="text-align:right;"> 0.476 </td>
   <td style="text-align:right;"> 1.60e-06 </td>
   <td style="text-align:left;"> IL4R,OSM,CISH,SOCS2,CD80,CDKN1A,MYC,PIM1,SOCS3,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:left;"> Cytokine Signaling in Immune system </td>
   <td style="text-align:right;"> 685 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.016 </td>
   <td style="text-align:right;"> 0.524 </td>
   <td style="text-align:right;"> 4.60e-06 </td>
   <td style="text-align:left;"> IL4R,OSM,CISH,EGR1,SOCS2,CD80,CDKN1A,MYC,PIM1,SOCS3,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 19 </td>
   <td style="text-align:left;"> Interleukin-4 and 13 signaling </td>
   <td style="text-align:right;"> 111 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.054 </td>
   <td style="text-align:right;"> 0.286 </td>
   <td style="text-align:right;"> 1.68e-05 </td>
   <td style="text-align:left;"> IL4R,OSM,CDKN1A,MYC,PIM1,SOCS3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 17 </td>
   <td style="text-align:left;"> Growth hormone receptor signaling </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.167 </td>
   <td style="text-align:right;"> 0.190 </td>
   <td style="text-align:right;"> 3.57e-05 </td>
   <td style="text-align:left;"> CISH,SOCS2,SOCS3,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> response to insulin </td>
   <td style="text-align:right;"> 255 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.027 </td>
   <td style="text-align:right;"> 0.189 </td>
   <td style="text-align:right;"> 1.06e-03 </td>
   <td style="text-align:left;"> CISH,EGR1,SOCS2,MYC,SOCS3,IRS2,PPARA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:left;"> Immune System </td>
   <td style="text-align:right;"> 2010 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.006 </td>
   <td style="text-align:right;"> 0.619 </td>
   <td style="text-align:right;"> 5.50e-03 </td>
   <td style="text-align:left;"> IL4R,OSM,CISH,EGR1,SOCS2,CD80,CDKN1A,MYC,PIM1,MB21D1,HSPA6,SOCS3,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> response to oxygen-containing compound </td>
   <td style="text-align:right;"> 1555 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.008 </td>
   <td style="text-align:right;"> 0.351 </td>
   <td style="text-align:right;"> 7.17e-03 </td>
   <td style="text-align:left;"> CPEB4,CISH,KLF9,EGR1,SOCS2,CDKN1A,MYC,PIM1,AQP3,RARG,SOCS3,IRS2,PPARA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:left;"> TFAP2 (AP-2) family regulates transcription of cell cycle factors </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.400 </td>
   <td style="text-align:right;"> 0.095 </td>
   <td style="text-align:right;"> 1.11e-02 </td>
   <td style="text-align:left;"> CDKN1A,MYC </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:left;"> Interleukin-7 signaling </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.083 </td>
   <td style="text-align:right;"> 0.143 </td>
   <td style="text-align:right;"> 1.38e-02 </td>
   <td style="text-align:left;"> CISH,SOCS2,IRS2 </td>
  </tr>
</tbody>
</table>

```r
# kable(head(term.induced[,c('p.value', 'term.name','intersection')], 10))
```

#### Untreated (WT) vs EPO with HMGB1 (EH)

##### Plots {.tabset .tabset-fade .tabset-pills}
###### MA-plot

```r
res_WT_EH <- results(wald.test, contrast = c("condition", "WT", "EH"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_WT_EH_ashr <- lfcShrink(wald.test, contrast = c("condition", "WT", "EH"), 
    res = res_WT_EV, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_WT_EH, xlim = xlim, ylim = ylim, main = "normal")
plotMA(res_WT_EH_ashr, xlim = xlim, ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 28-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_WT_EH$padj, breaks = 20, col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 29-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_WT_EH <- res_WT_EH[order(res_WT_EH$padj), ]
write.csv(as.data.frame(resOrdered_WT_EH), file = "condition_WT_EH_alpha0.1_results.csv")
WT_EH_p1 <- enhanced_volcano_plots("WT vs. EH (padj=0.05,log2fc=0.5)", res_WT_EH, 
    0.05, fc_cutoff)
WT_EH_p2 <- enhanced_volcano_plots("WT vs. EH (padj=0.1,log2fc=0.5)", res_WT_EH, 
    p_cutoff, fc_cutoff)
grid.arrange(WT_EH_p1, WT_EH_p2, nrow = 2, ncol = 1)
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
    print(plotCounts_gg(i, dds = wald.test, intgroup = c("condition")))
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
induced.sign <- rownames(res_WT_EV.df)[res_WT_EH.df$log2FoldChange >= 0.5 & 
    res_WT_EH.df$padj < p_cutoff]
# head(induced.sign) names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), ]
    # term.induced$p.value
    na.omit(term.induced)
    if (nrow(term.induced) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.induced)
    }
    kable(term.induced[1:rowlim, c("term.name", "term.size", "query.size", "overlap.size", 
        "recall", "precision", "p.value", "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis wit gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: induced gene functional analysis wit gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> Sulfide oxidation to sulfate </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.333 </td>
   <td style="text-align:right;"> 0.118 </td>
   <td style="text-align:right;"> 0.00963 </td>
   <td style="text-align:left;"> TST,SLC25A10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> Degradation of cysteine and homocysteine </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.154 </td>
   <td style="text-align:right;"> 0.118 </td>
   <td style="text-align:right;"> 0.04970 </td>
   <td style="text-align:left;"> TST,SLC25A10 </td>
  </tr>
</tbody>
</table>

```r
repressed.sign <- rownames(res_WT_EH.df)[res_WT_EH.df$log2FoldChange <= -0.5 & 
    res_WT_EH.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), ]
    na.omit(term.repressed)
    if (nrow(term.repressed) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.repressed)
    }
    kable(term.repressed[1:rowlim, c("term.name", "term.size", "query.size", 
        "overlap.size", "recall", "precision", "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: repressed genes functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 25 </td>
   <td style="text-align:left;"> Cytokine Signaling in Immune system </td>
   <td style="text-align:right;"> 685 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.015 </td>
   <td style="text-align:right;"> 0.476 </td>
   <td style="text-align:right;"> 6.82e-05 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,RELB,CISH,EGR1,SOCS2,IRF1,MYC,PIM1,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 24 </td>
   <td style="text-align:left;"> Immune System </td>
   <td style="text-align:right;"> 2010 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.007 </td>
   <td style="text-align:right;"> 0.714 </td>
   <td style="text-align:right;"> 7.65e-05 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,RELB,CISH,TNFAIP3,EGR1,SOCS2,IRF1,MYC,PIM1,NFKBIE,MB21D1,SLCO4C1,JUN,UBOX5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:left;"> Epstein-Barr virus infection </td>
   <td style="text-align:right;"> 197 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.036 </td>
   <td style="text-align:right;"> 0.304 </td>
   <td style="text-align:right;"> 8.57e-05 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,RELB,TNFAIP3,MYC,NFKBIE,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:left;"> HTLV-I infection </td>
   <td style="text-align:right;"> 253 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.028 </td>
   <td style="text-align:right;"> 0.304 </td>
   <td style="text-align:right;"> 4.54e-04 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,RELB,EGR1,MYC,ETS2,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 27 </td>
   <td style="text-align:left;"> Factor: E2F; motif: GGCGSG; match class: 1 </td>
   <td style="text-align:right;"> 10205 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.003 </td>
   <td style="text-align:right;"> 0.917 </td>
   <td style="text-align:right;"> 5.29e-04 </td>
   <td style="text-align:left;"> MSMO1,VMP1,KLF6,ICAM1,NFKBIA,RELB,ARRDC3,CISH,TNFAIP3,EGR1,SOCS2,ARL4A,IRF1,A4GALT,TMEM160,NINJ1,BHLHE40,MYC,PIM1,NFKBIZ,NFKBIE,RAB11FIP1,ETS2,ANKRD33B,MB21D1,KCNK5,STARD5,SLCO4C1,JUN,DDN,H1FX,UBOX5,TSC22D2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:left;"> TNF-alpha/NF-kappa B signaling complex (CHUK, KPNA3, NFKB2, NFKBIB, REL, IKBKG, NFKB1, NFKBIE, RELB, NFKBIA, RELA, TNIP2) </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 0.300 </td>
   <td style="text-align:right;"> 9.73e-04 </td>
   <td style="text-align:left;"> NFKBIA,RELB,NFKBIE </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:left;"> transcription factor activity, RNA polymerase II proximal promoter sequence-specific DNA binding </td>
   <td style="text-align:right;"> 399 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.020 </td>
   <td style="text-align:right;"> 0.222 </td>
   <td style="text-align:right;"> 1.32e-03 </td>
   <td style="text-align:left;"> KLF6,EGR1,IRF1,BHLHE40,MYC,ETS2,JUN,DDN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> RNA polymerase II regulatory region sequence-specific DNA binding </td>
   <td style="text-align:right;"> 634 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.014 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 4.17e-03 </td>
   <td style="text-align:left;"> KLF6,RELB,EGR1,IRF1,BHLHE40,MYC,ETS2,JUN,DDN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> RNA polymerase II regulatory region DNA binding </td>
   <td style="text-align:right;"> 635 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.014 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 4.23e-03 </td>
   <td style="text-align:left;"> KLF6,RELB,EGR1,IRF1,BHLHE40,MYC,ETS2,JUN,DDN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> double-stranded DNA binding </td>
   <td style="text-align:right;"> 832 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.012 </td>
   <td style="text-align:right;"> 0.278 </td>
   <td style="text-align:right;"> 4.47e-03 </td>
   <td style="text-align:left;"> KLF6,RELB,EGR1,IRF1,BHLHE40,MYC,ETS2,MB21D1,JUN,DDN </td>
  </tr>
</tbody>
</table>

```r
# kable(head(term.induced[,c('p.value', 'term.name','intersection')], 10))
```

***

#### Untreated (WT) vs HMGB1 (H) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_WT_H <- results(wald.test, contrast = c("condition", "WT", "H"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_WT_H_ashr <- lfcShrink(wald.test, contrast = c("condition", "WT", "H"), 
    res = res_WT_H, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_WT_H, xlim = xlim, ylim = ylim, main = "normal")
plotMA(res_WT_H_ashr, xlim = xlim, ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 34-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_WT_H$padj, breaks = 20, col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 35-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_WT_H <- res_WT_H[order(res_WT_H$padj), ]
write.csv(as.data.frame(resOrdered_WT_H), file = "condition_WT_H_alpha0.1_results.csv")
WT_H_p1 <- enhanced_volcano_plots("WT vs. H (padj=0.05,log2fc=0.5)", res_WT_H, 
    0.05, fc_cutoff)
WT_H_p2 <- enhanced_volcano_plots("WT vs. H (padj=0.1,log2fc=0.5)", res_WT_H, 
    p_cutoff, fc_cutoff)
grid.arrange(WT_H_p1, WT_H_p2, nrow = 2, ncol = 1)
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
    print(plotCounts_gg(i, dds = wald.test, intgroup = c("condition")))
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
induced.sign <- rownames(res_WT_H.df)[res_WT_H.df$log2FoldChange >= 0.5 & res_WT_H.df$padj < 
    p_cutoff]
# head(induced.sign) names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), ]
    na.omit(term.induced)
    if (nrow(term.induced) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.induced)
    }
    # term.induced$p.value
    if (length(term.induced$term.name) == 0) {
        cat("No terms found that have induced expression using the genelist")
    } else {
        kable(term.induced[1:rowlim, c("term.name", "term.size", "query.size", 
            "overlap.size", "recall", "precision", "p.value", "intersection")], 
            format.args = c(engeneer = TRUE, digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ") %>% 
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
                "responsive"))
    }
}
```

```
## No terms found that have induced expression using the genelist
```

```r
repressed.sign <- rownames(res_WT_H.df)[res_WT_H.df$log2FoldChange <= -0.5 & 
    res_WT_H.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), ]
    na.omit(term.repressed)
    if (nrow(term.repressed) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.repressed)
    }
    kable(term.repressed[1:rowlim, c("term.name", "term.size", "query.size", 
        "overlap.size", "recall", "precision", "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: repressed genes functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:left;"> Epstein-Barr virus infection </td>
   <td style="text-align:right;"> 197 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0.417 </td>
   <td style="text-align:right;"> 0.000472 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,TNFAIP3,NFKBIE,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> TNF signaling pathway </td>
   <td style="text-align:right;"> 108 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.037 </td>
   <td style="text-align:right;"> 0.333 </td>
   <td style="text-align:right;"> 0.001050 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,TNFAIP3,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> TNF-alpha/NF-kappa B signaling complex (CHUK, KPNA3, NFKB2, NFKBIB, REL, IKBKG, NFKB1, NFKBIE, RELB, NFKBIA, RELA, TNIP2) </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.167 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.003060 </td>
   <td style="text-align:left;"> NFKBIA,NFKBIE </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> CHUK-NFKB2-REL-IKBKG-SPAG9-NFKB1-NFKBIE-COPB2-TNIP1-NFKBIA-RELA-TNIP2 complex </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.167 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.003060 </td>
   <td style="text-align:left;"> NFKBIA,NFKBIE </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> regulation of DNA binding transcription factor activity </td>
   <td style="text-align:right;"> 400 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.015 </td>
   <td style="text-align:right;"> 0.300 </td>
   <td style="text-align:right;"> 0.005180 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,TNFAIP3,BHLHE40,NFKBIE,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> B cell receptor signaling pathway </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.043 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 0.008990 </td>
   <td style="text-align:left;"> NFKBIA,NFKBIE,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> nucleotide-binding oligomerization domain containing 1 signaling pathway </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.100 </td>
   <td style="text-align:right;"> 0.013800 </td>
   <td style="text-align:left;"> NFKBIA,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> Th1 and Th2 cell differentiation </td>
   <td style="text-align:right;"> 90 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.033 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 0.018900 </td>
   <td style="text-align:left;"> NFKBIA,NFKBIE,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> IL-17 signaling pathway </td>
   <td style="text-align:right;"> 92 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.033 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 0.020200 </td>
   <td style="text-align:left;"> NFKBIA,TNFAIP3,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:left;"> NF-kappa B signaling pathway </td>
   <td style="text-align:right;"> 93 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.032 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 0.020900 </td>
   <td style="text-align:left;"> ICAM1,NFKBIA,TNFAIP3 </td>
  </tr>
</tbody>
</table>

```r
# kable(head(term.induced[,c('p.value', 'term.name','intersection')], 10))
```

#### EPO with Vehicle (EV) vs EPO with HMGB1 (EH) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_EV_EH <- results(wald.test, contrast = c("condition", "EV", "EH"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_EV_EH_ashr <- lfcShrink(wald.test, contrast = c("condition", "EV", "EH"), 
    res = res_EV_EH, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_EV_EH, xlim = xlim, ylim = ylim, main = "normal")
plotMA(res_EV_EH_ashr, xlim = xlim, ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 40-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_EV_EH$padj, breaks = 20, col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 41-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_EV_EH <- res_EV_EH[order(res_EV_EH$padj), ]
write.csv(as.data.frame(resOrdered_EV_EH), file = "condition_EV_EH_alpha0.1_results.csv")
EV_EH_p1 <- enhanced_volcano_plots("EV vs. EH (padj=0.05,log2fc=0.5)", res_EV_EH, 
    0.05, fc_cutoff)
EV_EH_p2 <- enhanced_volcano_plots("EV vs. EH (padj=0.1,log2fc=0.5)", res_EV_EH, 
    p_cutoff, fc_cutoff)
grid.arrange(EV_EH_p1, EV_EH_p2, nrow = 2, ncol = 1)
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
    print(plotCounts_gg(i, dds = wald.test, intgroup = c("condition")))
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
induced.sign <- rownames(res_EV_EH.df)[res_EV_EH.df$log2FoldChange >= 0.5 & 
    res_EV_EH.df$padj < p_cutoff]
# head(induced.sign) names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), ]
    na.omit(term.induced)
    if (nrow(term.induced) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.induced)
    }
    # term.induced$p.value
    kable(term.induced[1:rowlim, c("term.name", "term.size", "query.size", "overlap.size", 
        "recall", "precision", "p.value", "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: induced gene functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> positive regulation of tau-protein kinase activity </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.087 </td>
   <td style="text-align:right;"> 0.0193 </td>
   <td style="text-align:left;"> EGR1,NAB2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> Attenuation phase </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.143 </td>
   <td style="text-align:right;"> 0.154 </td>
   <td style="text-align:right;"> 0.0334 </td>
   <td style="text-align:left;"> HSPA2,HSPA1A </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> regulation of tau-protein kinase activity </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.333 </td>
   <td style="text-align:right;"> 0.087 </td>
   <td style="text-align:right;"> 0.0481 </td>
   <td style="text-align:left;"> EGR1,NAB2 </td>
  </tr>
</tbody>
</table>

```r
repressed.sign <- rownames(res_EV_EH.df)[res_EV_EH.df$log2FoldChange <= -0.5 & 
    res_EV_EH.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), ]
    na.omit(term.repressed)
    if (nrow(term.repressed) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.repressed)
    }
    kable(term.repressed[1:rowlim, c("term.name", "term.size", "query.size", 
        "overlap.size", "recall", "precision", "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"), full_width = FALSE)
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>**Table: repressed genes functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 29 </td>
   <td style="text-align:left;"> TNF signaling pathway </td>
   <td style="text-align:right;"> 108 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.065 </td>
   <td style="text-align:right;"> 0.467 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,BCL3,ICAM1,NFKBIA,TNFAIP3,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25 </td>
   <td style="text-align:left;"> NF-kappa B signaling pathway </td>
   <td style="text-align:right;"> 93 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.065 </td>
   <td style="text-align:right;"> 0.400 </td>
   <td style="text-align:right;"> 9.00e-07 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,ICAM1,NFKBIA,RELB,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 32 </td>
   <td style="text-align:left;"> Epstein-Barr virus infection </td>
   <td style="text-align:right;"> 197 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.036 </td>
   <td style="text-align:right;"> 0.467 </td>
   <td style="text-align:right;"> 2.70e-06 </td>
   <td style="text-align:left;"> TRAF1,ICAM1,NFKBIA,RELB,TNFAIP3,NFKBIE,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 21 </td>
   <td style="text-align:left;"> I-kappaB/NF-kappaB complex </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.429 </td>
   <td style="text-align:right;"> 0.150 </td>
   <td style="text-align:right;"> 7.49e-05 </td>
   <td style="text-align:left;"> BCL3,NFKBIA,RELB </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> regulation of DNA binding transcription factor activity </td>
   <td style="text-align:right;"> 400 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.018 </td>
   <td style="text-align:right;"> 0.350 </td>
   <td style="text-align:right;"> 2.14e-04 </td>
   <td style="text-align:left;"> TRAF1,ICAM1,NFKBIA,TNFAIP3,BHLHE40,NFKBIE,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 23 </td>
   <td style="text-align:left;"> TNF-alpha/NF-kappa B signaling complex (CHUK, KPNA3, NFKB2, NFKBIB, REL, IKBKG, NFKB1, NFKBIE, RELB, NFKBIA, RELA, TNIP2) </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 0.375 </td>
   <td style="text-align:right;"> 3.12e-04 </td>
   <td style="text-align:left;"> NFKBIA,RELB,NFKBIE </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 38 </td>
   <td style="text-align:left;"> TNFR1-induced NFkappaB signaling pathway </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.115 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 7.68e-04 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> cellular response to tumor necrosis factor </td>
   <td style="text-align:right;"> 290 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.021 </td>
   <td style="text-align:right;"> 0.300 </td>
   <td style="text-align:right;"> 7.94e-04 </td>
   <td style="text-align:left;"> BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> response to tumor necrosis factor </td>
   <td style="text-align:right;"> 307 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.020 </td>
   <td style="text-align:right;"> 0.300 </td>
   <td style="text-align:right;"> 1.11e-03 </td>
   <td style="text-align:left;"> BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 39 </td>
   <td style="text-align:left;"> Regulation of TNFR1 signaling </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.091 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 1.61e-03 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,TNFAIP3 </td>
  </tr>
</tbody>
</table>

```r
# kable(head(term.induced[,c('p.value', 'term.name','intersection')], 10))
```

#### EPO with Vehicle (EV) vs HMGB1 (H) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_EV_H <- results(wald.test, contrast = c("condition", "EV", "H"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_EV_H_ashr <- lfcShrink(wald.test, contrast = c("condition", "EV", "H"), 
    res = res_EV_H, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_EV_H, xlim = xlim, ylim = ylim, main = "normal")
plotMA(res_EV_H_ashr, xlim = xlim, ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 46-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_EV_H$padj, breaks = 20, col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 47-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_EV_H <- res_EV_H[order(res_EV_H$padj), ]
write.csv(as.data.frame(resOrdered_EV_H), file = "condition_EV_H_alpha0.1_results.csv")
EV_H_p1 <- enhanced_volcano_plots("EV vs. H (padj=0.05,log2fc=0.5)", res_EV_H, 
    0.05, fc_cutoff)
EV_H_p2 <- enhanced_volcano_plots("EV vs. H (padj=0.1,log2fc=0.5)", res_EV_H, 
    p_cutoff, fc_cutoff)
grid.arrange(EV_H_p1, EV_H_p2, nrow = 2, ncol = 1)
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
    print(plotCounts_gg(i, dds = wald.test, intgroup = c("condition")))
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
induced.sign <- rownames(res_EV_H.df)[res_EV_H.df$log2FoldChange >= 0.5 & res_EV_H.df$padj < 
    p_cutoff]
# head(induced.sign) names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), ]
    # term.induced$p.value
    na.omit(term.induced)
    if (nrow(term.induced) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.induced)
    }
    kable(term.induced[1:rowlim, c("term.name", "term.size", "query.size", "overlap.size", 
        "recall", "precision", "p.value", "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: induced gene functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 36 </td>
   <td style="text-align:left;"> Growth hormone receptor signaling </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.208 </td>
   <td style="text-align:right;"> 0.147 </td>
   <td style="text-align:right;"> 3.30e-06 </td>
   <td style="text-align:left;"> CISH,SOCS2,SOCS3,SOCS1,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 32 </td>
   <td style="text-align:left;"> Cytokine Signaling in Immune system </td>
   <td style="text-align:right;"> 685 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.018 </td>
   <td style="text-align:right;"> 0.353 </td>
   <td style="text-align:right;"> 2.39e-04 </td>
   <td style="text-align:left;"> OSM,CISH,EGR1,SOCS2,CD80,CDKN1A,PIM1,EIF4A1,FOS,SOCS3,SOCS1,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 33 </td>
   <td style="text-align:left;"> Signaling by Interleukins </td>
   <td style="text-align:right;"> 463 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.022 </td>
   <td style="text-align:right;"> 0.294 </td>
   <td style="text-align:right;"> 3.92e-04 </td>
   <td style="text-align:left;"> OSM,CISH,SOCS2,CD80,CDKN1A,PIM1,FOS,SOCS3,SOCS1,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 34 </td>
   <td style="text-align:left;"> Interleukin-4 and 13 signaling </td>
   <td style="text-align:right;"> 111 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.054 </td>
   <td style="text-align:right;"> 0.176 </td>
   <td style="text-align:right;"> 4.01e-04 </td>
   <td style="text-align:left;"> OSM,CDKN1A,PIM1,FOS,SOCS3,SOCS1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30 </td>
   <td style="text-align:left;"> Jak-STAT signaling pathway </td>
   <td style="text-align:right;"> 162 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.043 </td>
   <td style="text-align:right;"> 0.206 </td>
   <td style="text-align:right;"> 4.96e-04 </td>
   <td style="text-align:left;"> OSM,CISH,SOCS2,CDKN1A,PIM1,SOCS3,SOCS1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:left;"> negative regulation of kinase activity </td>
   <td style="text-align:right;"> 272 </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.033 </td>
   <td style="text-align:right;"> 0.136 </td>
   <td style="text-align:right;"> 1.01e-03 </td>
   <td style="text-align:left;"> GADD45B,CISH,SOCS2,CDKN1A,DUSP2,SH3BP5L,SOCS3,SOCS1,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 31 </td>
   <td style="text-align:left;"> Prolactin signaling pathway </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.071 </td>
   <td style="text-align:right;"> 0.147 </td>
   <td style="text-align:right;"> 1.04e-03 </td>
   <td style="text-align:left;"> CISH,SOCS2,FOS,SOCS3,SOCS1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 37 </td>
   <td style="text-align:left;"> Factor: E2F; motif: GGCGSG; match class: 1 </td>
   <td style="text-align:right;"> 10205 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 0.005 </td>
   <td style="text-align:right;"> 0.815 </td>
   <td style="text-align:right;"> 1.05e-03 </td>
   <td style="text-align:left;"> KLF6,GADD45B,OSM,EIF5,EEF2K,RASD1,CISH,MREG,KLF9,EGR1,SOCS2,ARL4A,IKZF4,CDKN1A,SOX21,ID1,HSPA2,CHAC1,PIM1,ARRDC4,SIK1,RGS16,RHOB,FAM83A,KLF10,RAB11FIP1,DUSP2,EIF4A1,AGGF1,MB21D1,KCNK5,C10ORF10,NAB2,FAM102A,RHEBL1,MAT2A,B3GALNT1,FOS,KLF11,STARD5,HSPA6,SH3BP5L,CDK5R1,SOCS3,OSBP2,SOCS1,IRS2,HPDL,TSC22D2,ZNF470,NRARP,CCDC71L,PIGW </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35 </td>
   <td style="text-align:left;"> Interleukin-7 signaling </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.111 </td>
   <td style="text-align:right;"> 0.118 </td>
   <td style="text-align:right;"> 1.57e-03 </td>
   <td style="text-align:left;"> CISH,SOCS2,SOCS1,IRS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 17 </td>
   <td style="text-align:left;"> negative regulation of transferase activity </td>
   <td style="text-align:right;"> 298 </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.030 </td>
   <td style="text-align:right;"> 0.136 </td>
   <td style="text-align:right;"> 2.18e-03 </td>
   <td style="text-align:left;"> GADD45B,CISH,SOCS2,CDKN1A,DUSP2,SH3BP5L,SOCS3,SOCS1,IRS2 </td>
  </tr>
</tbody>
</table>

```r
repressed.sign <- rownames(res_EV_H.df)[res_EV_H.df$log2FoldChange <= -0.5 & 
    res_EV_H.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), ]
    na.omit(term.repressed)
    if (nrow(term.repressed) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.repressed)
    }
    kable(term.repressed[1:rowlim, c("term.name", "term.size", "query.size", 
        "overlap.size", "recall", "precision", "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: repressed genes functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:left;"> TNF signaling pathway </td>
   <td style="text-align:right;"> 108 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.065 </td>
   <td style="text-align:right;"> 0.412 </td>
   <td style="text-align:right;"> 1.00e-07 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,BCL3,ICAM1,NFKBIA,TNFAIP3,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 17 </td>
   <td style="text-align:left;"> NF-kappa B signaling pathway </td>
   <td style="text-align:right;"> 93 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.054 </td>
   <td style="text-align:right;"> 0.294 </td>
   <td style="text-align:right;"> 9.21e-05 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,ICAM1,NFKBIA,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> cellular response to tumor necrosis factor </td>
   <td style="text-align:right;"> 290 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.024 </td>
   <td style="text-align:right;"> 0.292 </td>
   <td style="text-align:right;"> 9.55e-05 </td>
   <td style="text-align:left;"> BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3,ZFP36L2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> response to tumor necrosis factor </td>
   <td style="text-align:right;"> 307 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.023 </td>
   <td style="text-align:right;"> 0.292 </td>
   <td style="text-align:right;"> 1.41e-04 </td>
   <td style="text-align:left;"> BIRC3,TNFRSF9,TRAF1,ICAM1,NFKBIA,TNFAIP3,ZFP36L2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> liver development </td>
   <td style="text-align:right;"> 139 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.036 </td>
   <td style="text-align:right;"> 0.208 </td>
   <td style="text-align:right;"> 1.47e-03 </td>
   <td style="text-align:left;"> NFKBIA,TNFAIP3,SLCO2B1,CITED2,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:left;"> TNFR1-induced NFkappaB signaling pathway </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.115 </td>
   <td style="text-align:right;"> 0.200 </td>
   <td style="text-align:right;"> 1.53e-03 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> hepaticobiliary system development </td>
   <td style="text-align:right;"> 142 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.035 </td>
   <td style="text-align:right;"> 0.208 </td>
   <td style="text-align:right;"> 1.64e-03 </td>
   <td style="text-align:left;"> NFKBIA,TNFAIP3,SLCO2B1,CITED2,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 21 </td>
   <td style="text-align:left;"> Regulation of TNFR1 signaling </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.091 </td>
   <td style="text-align:right;"> 0.200 </td>
   <td style="text-align:right;"> 3.20e-03 </td>
   <td style="text-align:left;"> BIRC3,TRAF1,TNFAIP3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:left;"> Epstein-Barr virus infection </td>
   <td style="text-align:right;"> 197 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0.294 </td>
   <td style="text-align:right;"> 3.60e-03 </td>
   <td style="text-align:left;"> TRAF1,ICAM1,NFKBIA,TNFAIP3,JUN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> tumor necrosis factor-mediated signaling pathway </td>
   <td style="text-align:right;"> 174 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.029 </td>
   <td style="text-align:right;"> 0.208 </td>
   <td style="text-align:right;"> 4.46e-03 </td>
   <td style="text-align:left;"> BIRC3,TNFRSF9,TRAF1,NFKBIA,TNFAIP3 </td>
  </tr>
</tbody>
</table>

```r
# kable(head(term.induced[,c('p.value', 'term.name','intersection')], 10))
```

#### EPO with HMGB1 (EH) vs HMGB1 (H) 

##### Plots {.tabset .tabset-fade .tabset-pills}

###### MA-plot

```r
res_EH_H <- results(wald.test, contrast = c("condition", "EH", "H"), alpha = p_cutoff, 
    pAdjustMethod = "BH")
res_EH_H_ashr <- lfcShrink(wald.test, contrast = c("condition", "EH", "H"), 
    res = res_EH_H, type = "ashr")
par(mfrow = c(1, 2))
plotMA(res_EH_H, xlim = xlim, ylim = ylim, main = "normal")
plotMA(res_EH_H_ashr, xlim = xlim, ylim = ylim, main = "ashr")
```

![](deseq2_analysis_files/figure-html/code block 52-1.png)<!-- -->

```r
invisible(dev.off())
```

###### Historgram of adjusted p-values

```r
hist(res_EH_H$padj, breaks = 20, col = "grey", main = "DESeq2 p-value distribution", 
    xlab = "DESeq2 P-value", ylab = "Number of genes")
```

![](deseq2_analysis_files/figure-html/code block 53-1.png)<!-- -->

###### Volacano-plot

```r
resOrdered_EH_H <- res_EH_H[order(res_EH_H$padj), ]
write.csv(as.data.frame(resOrdered_EH_H), file = "condition_EH_H_alpha0.1_results.csv")
EH_H_p1 <- enhanced_volcano_plots("EH vs. H (padj=0.05,log2fc=0.5)", res_EH_H, 
    0.05, fc_cutoff)
EH_H_p2 <- enhanced_volcano_plots("EH vs. H (padj=0.1,log2fc=0.5)", res_EH_H, 
    p_cutoff, fc_cutoff)
grid.arrange(EH_H_p1, EH_H_p2, nrow = 2, ncol = 1)
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
    print(plotCounts_gg(i, dds = wald.test, intgroup = c("condition")))
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
induced.sign <- rownames(res_EH_H.df)[res_EH_H.df$log2FoldChange >= 0.5 & res_EH_H.df$padj < 
    p_cutoff]
# head(induced.sign) names(term.induced)
if (identical(induced.sign, character(0))) {
    cat("No genes found that have induced expression")
} else {
    term.induced <- gprofiler(query = induced.sign, organism = "hsapiens")
    term.induced <- term.induced[order(term.induced$p.value), ]
    # term.induced$p.value
    na.omit(term.induced)
    if (nrow(term.induced) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.induced)
    }
    kable(term.induced[1:rowlim, c("term.name", "term.size", "query.size", "overlap.size", 
        "recall", "precision", "p.value", "intersection")], format.args = c(engeneer = TRUE, 
        digits = 3), caption = "**Table: induced gene functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>**Table: induced gene functional analysis with gProfileR. ** </caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> term.name </th>
   <th style="text-align:right;"> term.size </th>
   <th style="text-align:right;"> query.size </th>
   <th style="text-align:right;"> overlap.size </th>
   <th style="text-align:right;"> recall </th>
   <th style="text-align:right;"> precision </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> intersection </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> Jak-STAT signaling pathway </td>
   <td style="text-align:right;"> 162 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0.667 </td>
   <td style="text-align:right;"> 0.000146 </td>
   <td style="text-align:left;"> OSM,CISH,SOCS2,PIM1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> Cytokine Signaling in Immune system </td>
   <td style="text-align:right;"> 685 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.007 </td>
   <td style="text-align:right;"> 0.714 </td>
   <td style="text-align:right;"> 0.004300 </td>
   <td style="text-align:left;"> OSM,CISH,EGR1,SOCS2,PIM1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> Growth hormone receptor signaling </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.083 </td>
   <td style="text-align:right;"> 0.286 </td>
   <td style="text-align:right;"> 0.020700 </td>
   <td style="text-align:left;"> CISH,SOCS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> Signaling by Interleukins </td>
   <td style="text-align:right;"> 463 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.009 </td>
   <td style="text-align:right;"> 0.571 </td>
   <td style="text-align:right;"> 0.023100 </td>
   <td style="text-align:left;"> OSM,CISH,SOCS2,PIM1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> Interleukin-7 signaling </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.056 </td>
   <td style="text-align:right;"> 0.286 </td>
   <td style="text-align:right;"> 0.047100 </td>
   <td style="text-align:left;"> CISH,SOCS2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> RCP-Rab11 complex </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.049800 </td>
   <td style="text-align:left;"> RAB11FIP1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> EGR-EP300 complex </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.500 </td>
   <td style="text-align:right;"> 0.049800 </td>
   <td style="text-align:left;"> EGR1 </td>
  </tr>
</tbody>
</table>

```r
repressed.sign <- rownames(res_EH_H.df)[res_EH_H.df$log2FoldChange <= -0.5 & 
    res_EH_H.df$padj < p_cutoff]
if (identical(repressed.sign, character(0))) {
    cat("No genes found that have repressed expression")
} else {
    term.repressed <- gprofiler(query = repressed.sign, organism = "hsapiens")
    term.repressed <- term.repressed[order(term.repressed$p.value), ]
    na.omit(term.repressed)
    if (nrow(term.repressed) >= 10) {
        rowlim = 10
    } else {
        rowlim = nrow(term.repressed)
    }
    kable(term.repressed[1:rowlim, c("term.name", "term.size", "query.size", 
        "overlap.size", "recall", "precision", "p.value", "intersection")], 
        format.args = c(engeneer = TRUE, digits = 3), caption = "**Table: repressed genes functional analysis with gProfileR. ** ") %>% 
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", 
            "responsive"))
}
```

```
## No genes found that have repressed expression
```

```r
# kable(head(term.induced[,c('p.value', 'term.name','intersection')], 10))
```

Differential Gene Expression (DGE) Analysis 
================
Tocora et al., 2026

The following pipeline is modified from "Analyzing RNA-seq data with DESeq2" by Love et al (2022) available here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

### Installing and loading packages
```{r Install and Load Packages}
### DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")

### EnhacedVolcano (https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) and (https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)
devtools::install_github('kevinblighe/EnhancedVolcano')
#OR
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
  BiocManager::install('EnhancedVolcano')
  
### Pheatmap (https://cran.r-project.org/web/packages/pheatmap/index.html)
install.packages("pheatmap")

### tidyverse (https://www.tidyverse.org/)
install.packages("tidyverse")

### RColorBrewer (https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf)
install.packages("RColorBrewer")

### UpSetR (https://github.com/hms-dbmi/UpSetR)
install.packages("UpSetR")

### ComplexHeatmap (https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library("DESeq2")
library(tidyverse)
library(RColorBrewer)
library("pheatmap")
library(lubridate)
library(cowplot)
library(lme4)
library(png)
library(car)
library(magick)
library(lm.beta)
library(UpSetR)
library(ComplexHeatmap)
````

### Imput data for activity comparison (High- vs Low-activity Bodyguards) 
We manually obtained a matrix of read counts from the files obtained in Salmon and imported this matrix under "countdata". I also imported the sample information table "coldata". The "countdata" matrix has the samples information in columns organized in the same order as the rows in the "coldata" table.  then, construct a DESeqDataSet with the DESeqDataSetFromMatrix() function.  

```{r Install and Load Packages}
countdata <- read.csv("activity_count_matrix.csv", header = TRUE, row.names = 1)
coldata <- read.delim("activity_samples.txt", header = TRUE, row.names = 1)
all(rownames(coldata) == colnames(countdata))
### Need TRUE as output 
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~activity.level)
ddsFullCountTable
```

The DESeqDataSet got as output is described as follows: 

- class: DESeqDataSet 
- dim: 67613 14 
- metadata(1): version
- assays(1): counts
- rownames(67613): g0 g1 ... g9998 g9999
- rowData names(0):
- colnames(14): sample1 sample2 ... sample17 sample19
- colData names(5): plant.id code treatment attack.avg activity.level

### Differential expression analysis for the activity comparison
Differential gene expression analysis based on the negative binomial distribution described in Love et al (2014). For mor information check the Methods section of the paper. 

```{r Install and Load Packages}
dds <- ddsFullCountTable
dds$activity.level <- relevel(dds$activity.level, "low")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
````
The output obtained "res" containes: 
- log2 fold change (MLE): activity.level high vs low 
- Wald test p-value: activity.level high vs low 
- DataFrame with 67613 rows and 6 columns

### Extracting significant differentially expressed genes in the activity comparison
The following pipeline modifies code from https://hbctraining.github.io/Intro-to-R-with-DGE/lessons/B1_DGE_visualizing_results.html
Please consider that thet The lfc.cutoff is set to 0.58, so since we are working with log2 fold changes so this translates to an actual fold change of 1.5; and the padj.cutoff is set to 0.01. The final databases "Activity_sigOE.csv" and " are available in the folder "Data"  

```{r Install and Load Packages}
res_ordered = res[order(res$pvalue),]
head(res_ordered)
summary(res_ordered)
sum(res_ordered$padj < 0.01, na.rm=TRUE)
### Setting up parameters: 
padj.cutoff <- 0.01
lfc.cutoff <- 0.58
threshold <- res_ordered$padj < padj.cutoff & abs(res_ordered$log2FoldChange) > lfc.cutoff ### vector containing genes that meet our criteria (parameters stated above). 
length(which(threshold)) ### the vector has a length equal to the total number of significant genes in the dataset.
res_ordered$threshold <- threshold 
sigOE <- data.frame(subset(res_ordered, threshold==TRUE)) ## subsetting significant genes
write.csv(sigOE, file="Activity_sigOE.csv")
normalized_counts <- counts(dds, normalized=T) ###normalized counts
norm_OEsig <- normalized_counts[rownames(sigOE),] ###Extract normalized expression for significant genes
write.csv(norm_OEsig, file="Activity_norm_OEsig.csv")
```

### _Figure 2: Plotting Significant genes_
In this section, I extract the normalized values of all the significant genes and plot a heatmap of their expression using pheatmap()

### Figure 2A. 
```{r Install and Load Packages}
### Run pheatmap
heat.colors <- brewer.pal(6, "Greys") ### Set a color palette

annotation <- data.frame(activity.level=coldata[,'activity.level'], 
                     row.names=rownames(coldata))
my_colour = list(activity.level = c(high = "#BA55D3", low = "#FFD700"))

p1 <- pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,annotation= annotation, border_color=NA, fontsize = 10, scale="row", fontsize_row = 10, height=20, annotation_colors = my_colour, legend_labels = F, annotation_names_row = F, annotation_names_col = F)
```

### Figure 2B-C
```{r echo=TRUE}
res2 <- results(dds, tidy=TRUE, contrast=c("activity.level", "high", "low")) %>%
  arrange(padj, pvalue) %>%
  tibble::as_tibble()
res2

# Figure 2B.

#Metalloproteinase (g16223)
##Where is the gene we want to plot? 
which(grepl("g16223", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[45]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
p2 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="Bodyguard Activity", y="Expression (log normalized counts)", title = "Metalloproteinase", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_fill_manual(values = alpha(c( "#ff1493", "#0000ff"), 0.4)) +  geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="padj = 0.006", color="black")

#check ylim in the tcounts dataframe
###For no legend: theme(legend.position="none")
##Where is the significant value? 
sigOE["g16223",]

#Figure 2C. 
#Chymotrypsin (g7229)
##Where is the gene we want to plot? 
which(grepl("g7229", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[23]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
p3 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="Bodyguard Activity", y="Expression (log normalized counts)", title = "Chymotrypsin 1-like", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +scale_fill_manual(values = alpha(c( "#ff1493", "#0000ff"), 0.4)) +  geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=13, label="padj = 0.001", color="black")
##Where is the significant value? 
sigOE["g7229",]
```

### Caste_analysis (all workers)
With the following code I perform the analysis for the caste comparison. Please consider that I used as imput files the "caste_count_matrix.csv", and the samples information table "caste_samples.csv". The final databases "Caste_sigOE.csv" and "Caste_norm_OEsig.csv" are in the data folder.  

```{r Install and Load Packages}
countdata <- read.csv("CASTEgene_count_matrix_collapsedTranscriptomeWithDecoys.csv", header = TRUE, row.names = 1)
coldata <- read.delim("caste_modified.txt", header = TRUE, row.names = 1)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~treatment)
ddsFullCountTable
dds <- ddsFullCountTable
dds$treatment <- relevel(dds$treatment, "brood")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
```

ddsFullCountTable: 
- class: DESeqDataSet 
- dim: 67613 12 
- metadata(1): version
- assays(1): counts
- rownames(67613): g0 g1 ... g9998 g9999
- rowData names(0):
- colnames(12): sample4 sample5 ... sample19 sample20
- colData names(13): plant.id treatment.code ... date time

res: 
- log2 fold change (MLE): treatment guard vs brood 
- Wald test p-value: treatment guard vs brood 
- DataFrame with 67613 rows and 6 columns

```{r Install and Load Packages}
### Extracting significant differentially expressed genes
res_ordered = res[order(res$pvalue),]
head(res_ordered)
summary(res_ordered)
sum(res_ordered$padj < 0.01, na.rm=TRUE)
### Setting up parameters: 
padj.cutoff <- 0.01
lfc.cutoff <- 0.58
threshold <- res_ordered$padj < padj.cutoff & abs(res_ordered$log2FoldChange) > lfc.cutoff ### vector containing genes that meet our criteria (parameters stated above). 
length(which(threshold)) ### the vector has a length equal to the total number of significant genes in the dataset.
res_ordered$threshold <- threshold 
sigOE <- data.frame(subset(res_ordered, threshold==TRUE)) ## subsetting significant genes
write.csv(sigOE, file="Caste_sigOE.csv")
normalized_counts <- counts(dds, normalized=T) ###normalized counts
norm_OEsig <- normalized_counts[rownames(sigOE),] ###Extract normalized expression for significant genes
write.csv(norm_OEsig, file="Caste_norm_OEsig.csv")
```

### Figure 2D - UpSet Plot

```{r echo=TRUE}
setwd("~/Documents/Trancriptomics A. octoarticulatus/Graphs")
###install.packages("ComplexUpset")
library(ggplot2)
library(ComplexUpset)
#mydata = read.csv("DGE_geneSummary.csv")
mydata = read.csv("DGE_Summary.csv")
Analysis = colnames(mydata)[9:11]


p4 <- upset(
  mydata,
  Analysis,
  width_ratio = 0.15,
  min_size = 5,
  
  # Drop legacy theme warnings
  themes = list(default = theme()),
  
  # 1. TOP PANEL: Intersection sizes (Black bars, white background)
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = 'black_bars')
    ) + 
    scale_fill_manual(values = c('black_bars' = '#000000'), guide = 'none') +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "#e0e0e0", linetype = "solid"),
      panel.grid.minor = element_blank()
    )
  ),
  
  # 2. SIDE PANEL: Set sizes (Black bars, white background)
  set_sizes = upset_set_size(
    mapping = aes(fill = 'black_bars')
  ) + 
  scale_fill_manual(values = c('black_bars' = '#000000'), guide = 'none') +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "#e0e0e0", linetype = "solid"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ),
  
  # 3. BOTTOM PANEL: Matrix Grid
  matrix = intersection_matrix() + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(size = 11, face = "bold", color = "#000000")
    )
) & 
# Apply global aesthetic white-out adjustments across the combined framework
theme(
  plot.background = element_rect(fill = "white", color = NA),
  text = element_text(family = "sans")
)


# Calculate the total size for each group/set
set_sizes <- colSums(mydata[, Analysis, drop = FALSE])

# Print the result
print(set_sizes)
```

### Supplementary: Visualization of top significant genes for Activity and Caste comparisons
In this section I plotted the top differentially expressed genes considering a log2FoldChange equal to 2 and padj value of 10e-6 as cutoff parameters. I used the R package "EnhacedVolcano" (Blighe, et al., 2018). Please consider that the log transformed adjusted p-values are plotted on the y-axis and log2 fold change values on the x-axis. 

### Activity Comparison

```{r Install and Load Packages}
keyvals <- ifelse(
    res$log2FoldChange < -2 & res$padj < 10e-6, 'royalblue',
      ifelse(res$log2FoldChange > 2 & res$padj < 10e-6, 'gold',
        'black')) 
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'Upregulated high active bodyguards'
  names(keyvals)[keyvals == 'black'] <- 'NS/Only Log2 FC/Only p-value'
  names(keyvals)[keyvals == 'royalblue'] <- 'Upregulated low active bodyguards'
  
 ActivityVolcanoPlot <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res)[which(names(keyvals) %in% c('Upregulated high active bodyguards', 'Upregulated low active bodyguards'))],
    xlab = bquote(~Log[2]~ 'fold change'),
    title = 'Activity Analysis',
    pCutoff = 10e-6,
    FCcutoff = 2.0,
    pointSize = 4,
    labSize = 4.5,
    boxedLabels = TRUE,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'top',
    legendLabSize = 15,
    legendIconSize = 5.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    arrowheads = TRUE,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black')
ActivityVolcanoPlot
```
### Caste Comparison

```{r Install and Load Packages}
   keyvals <- ifelse(
    res$log2FoldChange < -2 & res$padj < 10e-6, 'royalblue',
      ifelse(res$log2FoldChange > 2 & res$padj < 10e-6, 'gold',
        'black')) 
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'Upregulated bodyguards'
  names(keyvals)[keyvals == 'black'] <- 'NS/Only Log2 FC/Only p-value'
  names(keyvals)[keyvals == 'royalblue'] <- 'Upregulated brood-care workers'
  
  CasteVolcanoPlot <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res)[which(names(keyvals) %in% c('Upregulated bodyguards', 'Upregulated brood-care workers'))],
    xlab = bquote(~Log[2]~ 'fold change'),
    title = 'Caste Analysis',
    pCutoff = 10e-6,
    FCcutoff = 2.0,
    xlim = c(-10,10),
    pointSize = 3.5,
    labSize = 3.5,
    labCol = 'black',
    boxedLabels = TRUE,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'top',
    legendLabSize = 15,
    legendIconSize = 5.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    arrowheads = TRUE,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black)
    
 CasteVolcanoPlot
```

## __References__

1. DESeq2: Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.
2. EnhacedVolcano: Blighe, K, S Rana, and M Lewis. 2018. “EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling.” https://github.com/kevinblighe/EnhancedVolcano.

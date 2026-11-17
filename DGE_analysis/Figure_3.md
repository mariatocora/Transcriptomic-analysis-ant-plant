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
### _Figure 3. Differences in viral transcript expression between high- and low-activity bodyguards across colonies_

```{r echo=TRUE}
res2 <- results(dds, tidy=TRUE, contrast=c("activity.level", "high", "low")) %>%
  arrange(padj, pvalue) %>%
  tibble::as_tibble()
res2

##Narnavirus (g9141) (AoV1)
which(grepl("g9141", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[7]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g9141<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-1", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 6.45E-11", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g9141

#Putative A. octoarticulatus virus 2 (g6982_i4) AoV2
##Where is the gene we want to plot? 
which(grepl("g6982", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[4]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g6982<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-2", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 7.8E-13", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g6982

#Putative A. octoarticulatus virus 3 (g2275_i0) (AoV-3)
##Where is the gene we want to plot? 
which(grepl("g2275", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[25]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g2275 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-3", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 0.001", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g2275
#check ylim in the tcounts dataframe
###For no legend: theme(legend.position="none")

#Putative A. octoarticulatus virus 4 (g5994_i0) AoV4
##Where is the gene we want to plot? 
which(grepl("g5994", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[17]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g5994<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-4", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 0.0001", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g5994


#Putative A. octoarticulatus virus 4 (g1894_i0) AoV-5
##Where is the gene we want to plot? 
which(grepl("g1894", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[11]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g1894<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-5", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 1.4E-6", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g1894

###Putative A. octoarticulatus virus 5 (g2855_i0) AoV6
which(grepl("g2855", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[15]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g2855<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-6", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 8.6E-5", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g2855

###Putative A. octoarticulatus virus 6 (g3684_i0) AoV-7
which(grepl("g3684", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[24]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g3684<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-7", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 0.001", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g3684

##Putative A. octoarticulatus virus 7 (g3814_i0) AoV-8
which(grepl("g3814", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[38]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g3814<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level", outliers = FALSE) +  labs(x="Activity level", y="Expression (log normalized counts)", title = "AoV-8", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "bold")) +  scale_fill_manual(values = alpha(c("#ff1493", "#0000ff"), 0.4)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="padj = 0.004", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("low" = "#ff1493", "high"  = "#0000ff"))
g3814
```

### Multipanel Virus A-H
```{r echo=TRUE}
up_row <- plot_grid(g9141, g6982, g2275, g5994, labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol = 4)
up_row
bottom_row <- plot_grid(g1894, g2855, g3684, g3814, labels = c('E', 'F', 'G','H'), label_size = 12, ncol = 4)
bottom_row
Virus_figures <- plot_grid(up_row, bottom_row, label_size = 12, ncol = 1)
Virus_figures
#Save plot
###save_plot("Figure6New.pdf", Figure_3, base_width=10, base_height=7)
```
### Pannel I and J Caste Comparison data
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

### Hymenoptaecin and Metallopeptidase

```{r echo=TRUE}
res2 <- results(dds, tidy=TRUE, contrast=c("activity.level", "high", "low")) %>%
  arrange(padj, pvalue) %>%
  tibble::as_tibble()
res2

###g68084_hymenoptaecin 
which(grepl("g68084", res2$row))

goi <- res2$row[12129]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g68084 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Hymenoptaecin", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 10)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=8, label="p = n.s.", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = alpha(c("#0000ff", "#ffd500", "#ff1493"), 0.4)) + annotate(geom="text", x =2, y=10, label="p = 2.67E-07", color="black") + annotate(geom="text", x =2.5, y=8, label="p = ns.s", color="black") + geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("bood" = "#ffd500", "low" = "#ff1493", "high"  = "#0000ff")) ##BvsH(p = 9.89E-07)##BvsL(ns)#HvsL(ns)
g68084


###g9196_metallopeptidase(3)
which(grepl("g9196", res2$row))

goi <- res2$row[75]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g9196 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Metallopeptidase", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5))+ theme(legend.position="none") + annotate(geom="text", x =1.5, y=8, label="p = 8.018E-06", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = alpha(c("#0000ff", "#ffd500", "#ff1493"), 0.4)) + annotate(geom="text", x =2, y=10, label="p = 2.675E-07", color="black") + annotate(geom="text", x =2.5, y=8, label="p = n.s.", color="black") +  geom_jitter(aes(color = activity.level), width = 0.2, size = 2) + scale_color_manual(values = c("bood" = "#ffd500", "low" = "#ff1493", "high"  = "#0000ff")) ##BvsH(p = 9.89E-07)##BvsL(ns)#HvsL(ns)
g9196

genes <- plot_grid(g68084, g9196, labels = c('I', 'J'), label_size = 12)
genes
```
### Final Panel

```{r echo=TRUE}
Figure_3 <- plot_grid(Virus_figures, genes, label_size = 12, ncol = 1, nrow = 2, rel_heights = c(2,1), rel_widths = c(1,1))

Figure_3

#Save plot
###save_plot("Figure6New.pdf", Figure_3, base_width=10, base_height=7)
```



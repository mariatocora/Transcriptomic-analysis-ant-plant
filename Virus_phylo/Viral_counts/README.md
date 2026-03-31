### Counts

```{r load packages, warning=FALSE, message=FALSE}
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggplot2)

##### RPKM
RPKM <- read.csv("RPKM_counts.csv")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Correlation <- read.csv("Correlation.csv")

p1 <- ggplot(data = Correlation, aes(x=Broodcare, y=Bodyguard, label = Colony)) + geom_point(aes(color=Virus, size=2)) + scale_color_manual(values = cbp2) + labs( x = "RPKM Broodcare", y = "RPKM Bodyguard") + theme_cowplot(12)

### genes
genes <- read.csv("bbsplit_genes_counts.csv")
p2 <- ggplot(data = genes, aes(x=gene, y=bbsplit, fill = Activity)) + geom_boxplot() + scale_fill_manual(values = c("gold", "#0000ff", "#ff1493")) + labs( x = "Gene", y = "% Unambiguous mapped reads") + theme_cowplot(12)

##virus
virus <- read.csv("bbsplit_virus_counts.csv")
virus$bbsplit <- as.numeric(virus$bbsplit)
p3 <- ggplot(data = virus, aes(x=Sample, y=bbsplit, color=Virus)) + 
    geom_point(aes(color=Virus, size=2)) + labs( x = "Sample", y = "% unambiguous mapped reads") + facet_grid(~Caste, scale = "free") + geom_hline(yintercept = 0.001, color = "red", linetype = "dotted") + scale_color_manual(values = cbp2) + theme_half_open(12) +  panel_border() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

upper_row <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
pannel <- plot_grid(upper_row, p3, labels = c('', 'C'), label_size = 12, ncol = 1)
```


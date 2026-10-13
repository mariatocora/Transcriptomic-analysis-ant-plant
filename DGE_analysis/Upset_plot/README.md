Upset Plot
==================================================
Maria Tocora 

The following code is modified from https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html#upset-plots-as-heatmaps. Tables with the number of genes per every differential expression analysis category are provided below under every comparison. The UpSet plots group co-occurring variables into sets and shows you a bar chart of their frequency (Lex et al., 2014; Conway, et al., 2017). The databases are available in the "Data" folder. 

## Install and load packages 
```{r Install and Load Packages}
### UpSetR (https://github.com/hms-dbmi/UpSetR)
install.packages("UpSetR")

### ComplexHeatmap (https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(UpSetR)
library(ComplexHeatmap)
```

### UpSet Plot
```
###UpSetPlot
```{r echo=TRUE}
setwd("~/Documents/Trancriptomics A. octoarticulatus/Graphs")
###install.packages("ComplexUpset")
library(ggplot2)
library(ComplexUpset)
mydata = read.csv("DGE_geneSummary.csv")
Analysis = colnames(mydata)[2:4]


upset(
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
```

## __References__

1. UpSetR: Conway, J. R., Lex, A., & Gehlenborg, N. (2017). UpSetR: an R package for the visualization of intersecting sets and their properties. Bioinformatics (Oxford, England), 33(18), 2938–2940. https://doi.org/10.1093/bioinformatics/btx364
2. UpSetR: Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister, H. (2014). UpSet: Visualization of Intersecting Sets. IEEE transactions on visualization and computer graphics, 20(12), 1983–1992. https://doi.org/10.1109/TVCG.2014.2346248
3. ComplexHeatmap: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.

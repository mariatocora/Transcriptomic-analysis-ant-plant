### Code to make Fig. S5
```{r setup, include=FALSE}
setwd("~/Documents/Trancriptomics A. octoarticulatus/AllData/Virus")
PoissonGamma <- read.csv("PoissonGamma.csv")
data <- PoissonGamma[, c("Species.1", "Species.2", "similarity")]
df1 <- as.data.frame.matrix(xtabs(similarity ~ ., data))
UpperT <- df1[upper.tri(df1, diag = TRUE)]
n <- dim(df1)[[1]]
Row <- unlist(lapply(seq_len(n), FUN = seq_len))
Column <- rep(seq_len(n), seq_len(n))
#### Match row and column names 
names <- c(rownames(df1))
numbers <- c(1:74)
Table <- data.frame(names, numbers)
#### Row
Df2 <- data.frame(UpperT, Row, Column)
x <- Df2$Row
Table$names[match(x, Table$numbers)]
Df2$Row <- Table$names[match(x, Table$numbers)]
### column
x <- Df2$Column
Table$names[match(x, Table$numbers)]
Df2$Column <- Table$names[match(x, Table$numbers)]
### Final
cleaned_df <- na.omit(Df2)

### Plot
PoissonGamma_plot <- ggplot(data = cleaned_df, aes(Row, Column, fill = UpperT))+ geom_tile(color = "white") + scale_fill_gradient2(low = "#5BA4BB", mid = "#FDFD96", high = "#FC979B", midpoint = 0.5, name="% AA similarity") + labs( x = "Species", y = "Species") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + coord_fixed() + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

###loadedNamespaces()
 [1] "methods"      "sass"         "graphics"     "generics"     "tidyr"       
 [6] "stringi"      "lattice"      "hms"          "digest"       "magrittr"    
[11] "RColorBrewer" "evaluate"     "grid"         "timechange"   "grDevices"   
[16] "fastmap"      "plyr"         "jsonlite"     "Matrix"       "tidyverse"   
[21] "mgcv"         "gaston"       "purrr"        "scales"       "textshaping" 
[26] "jquerylib"    "cli"          "crayon"       "rlang"        "cowplot"     
[31] "bit64"        "munsell"      "splines"      "cachem"       "withr"       
[36] "utils"        "yaml"         "stats"        "parallel"     "tools"       
[41] "base"         "reshape2"     "tzdb"         "dplyr"        "colorspace"  
[46] "ggplot2"      "forcats"      "vctrs"        "R6"           "lifecycle"   
[51] "lubridate"    "stringr"      "bit"          "vroom"        "ragg"        
[56] "pkgconfig"    "RcppParallel" "bslib"        "pillar"       "gtable"      
[61] "Rcpp"         "glue"         "systemfonts"  "xfun"         "tibble"      
[66] "tidyselect"   "rstudioapi"   "knitr"        "farver"       "datasets"    
[71] "htmltools"    "nlme"         "rmarkdown"    "labeling"     "readr"       
[76] "compiler"  
```

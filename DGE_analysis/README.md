Differential Gene Expression (DGE) Analysis 
================
Tocora et al., 
2026

For transcript quantification, we implemented the Salmon v1.9.0 pipeline (Patro et al., 2017). We built a Salmon index based on a decoy-aware transcriptome to avoid false mapping (Srivastava et al., 2020). We used the genome of a close relative of A. octoarticulatus, Wasmannia auropunctata (GenBank: GCA_000956235.1) as the ‘decoy.’ Then we aligned the RNA-Seq reads to the Salmon index. In this step, we implemented the Gibbs sampling option (--numGibbsSamples 20) and reduced autocorrelation in the Gibbs samples (thinningFactor 100). Finally, we obtained transcript-level counts for each sample with tximport (Soneson et al., 2015). 

To identify differentially expressed genes (DEGs), we used the DESeq2 workflow (Love et al., 2014) in R version 4.2.1 (R Core Team, 2020). We extracted the final list of significant genes with the results function and set a false discovery rate (FDR) < 0.01 and fold change > 1.5 (Log2FC value of 0.58). We made heatmaps with the R library pheatmap. 


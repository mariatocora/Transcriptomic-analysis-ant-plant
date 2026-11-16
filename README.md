Data analyses for Tocora et al: RNA-Seq finds viruses associated with mutualist quality in the Amazonian plant-ant Allomerus octoarticulatus
================
Tocora et al., 2026

### Abstract: 
Ant-plant symbioses are classic examples of mutualism in which ant “bodyguards” defend myrmecophytic plants against enemies in exchange for nest sites and often food. In these interactions, variation in ant defense has cascading effects on the fitness of their host plants, as plants that are less well protected by ants may experience higher herbivory. However, only a handful of studies have explored the molecular mechanisms behind how much ants invest in cooperation in their mutualisms with plants. Here, we used RNA-Seq to profile the transcriptomes of Allomerus octoarticulatus ant workers, which aggressively defend the Amazonian plant Cordia nodosa against herbivores, but to varying degrees. Field behavioral assays with herbivores in the Peruvian Amazon showed striking variation among colonies in the relative zeal with which A. octoarticulatus workers defend their host plant. Highly effective and ineffective bodyguards differed in their gene expression profiles, which revealed virus-like transcript abundance significantly associated with variation in ant bodyguarding behavior. Eight new positive-sense single-stranded RNA viruses were differentially expressed between colonies with high- or low-quality bodyguards. Gene expression, including that of immunity-related genes, also differed between broodcare workers and bodyguard ants, suggesting bodyguarding is a distinct worker task. Viruses may influence ant cooperation with plants in ant-plant mutualisms. 


### Description of pipeline, relevant folders in parentheses:
 1.  De novo transcriptome (./De_novo_transcriptome/) 
  1.1  Samples alignment to a Decoy-aware transcriptome and counts quantification. Based on Patro et al (2017)
 2.  De novo transcriptome annotation (./De_novo_transcriptome/Annotation)
 3.  Differential gene expression analysis (./DGE_analysis). Based on Love et al (2014) 
 4. Phylogenetic analysis of virus (./Virus_phylo)
 5. Field data analysis (.FieldData)
 6. References (./references): References for methods used in the analyses. 

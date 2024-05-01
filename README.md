# MiMedSurv

Title: MiMedSurv: A Unified Cloud Platform for Microbiome Causal Mediation Analysis with Survival Responses

Version: 1.0.0

Description: MiMedSurv (microbiome mediation analysis with survival responses) is a unified web cloud computing platform for comprehensive microbiome mediation analysis with survival (i.e., time-to-event) responses. MiMedSurv surveys the roles of the microbiome as a mediator between a treatment (e.g., medical intervention, environmental exposure) and a survival response on the host’s health or disease (e.g., time-to-disease, time-to-cure). The main features of MiMedSurv are as follows. First, MiMedSurv conducts basic exploratory non-mediational survival analysis, not involving microbiome, to survey the disparity in survival time between medical treatments (e.g., treatment vs. placebo, new treatment vs. old treatment) / environmental exposures (e.g., rural vs. urban, smoking vs. non-smoking). (see Non-Mediational Analysis). Second, MiMedSurv identifies the mediating roles of the microbiome in various aspects (see Mediational Analysis): (i) as a microbial ecosystem using ecological measures (e.g., alpha- and beta-diversity indices) (see Community-level Analysis) and (ii) as individual microbial taxa in various hierarchies (e.g., phyla, classes, orders, families, genera, species) (see Taxonomy-level Analysis). We also stress that MiMedSurv can conduct covariate-adjusted analysis to control for potential confounding factors (e.g., age, sex) to enhance the causality of the results especially for observational studies. MiMedSurv also provides user-friendly data preprocessing and analytic modules and makes nice visualizations.

NeedsCompilation: No

Depends: R(≥ 4.1.0)

Imports: Bioconductor ('BiocParallel', 'biomformat', 'phyloseq'); CRAN ('betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zip', 'bda', 'mediation', 'survival', 'survminer', 'coin'); GitHub ('DACT')

License: GPL 1, GPL 2 

**URLs**: Web Server (http://mimedsurv.micloud.kr), GitHub (http://github.com/yj7599/mimedsurvgit) 

**Maintainer**: Hyojung Jang (hyojung.jang@northwestern.edu)

## Reference

* Jang, H., Koh, H. A unified web cloud computing platform MiMedSurv for microbiome causal mediation analysis with survival responses (In Review)

## Descriptions on Data Processing & Analytic Modules

**(1) Data Processing: Data Input** - Users first need to upload their microbiome data with (1) four components: a feature (operational taxonomic unit (OTU) or amplicon sequence variants (ASV)) table, a taxonomic table, a phylogenetic tree, and metadata or (2) three components: a feature (OTU or ASV) table, a taxonomic table, and metadata. If users upload their microbiome dataset with three components without a phylogenetic tree, only non-phylogenetic community-level (alpha- and beta-diversity) analyses will later be performed. Users can start with downloading the example microbiome data from the Example Data section. The data are stored in a widely used unified format, called phyloseq, that can efficiently combine all essential microbiome data components. Alternatively, users can also employ them all individually. In this module, we also described all the instructions with example codes to check compatible data formats so that users can prepare microbiome data easily. 

**(2) Data Processing: Quality Control** - Users can perform quality controls with respect to (i) a microbial kingdom of interest (default: Bacteria), (ii) a minimum library size (i.e., total read count) for the study subjects to be retained (default: 3,000), (3) a minimum mean proportion for the features (OTUs or ASVs) to be retained (default: 0.02%), and (4) errors in taxonomic names to be removed. 

**(3) Data Processing: Data Transformation** - For the community-level analysis, users can compute nine alpha diversity indices (i.e., Observed, Shannon, Simpson, Inverse Simpson, Fisher, Chao1, abundance-based coverage estimator (ACE), incidence-based coverage estimator (ICE), phylogenetic diversity (PD)) and five beta diversity indices (i.e., Jaccard dissimilarity, Bray-Curtis dissimilarity, Unweighted UniFrac distance, Generalized UniFrac distance, Weighted UniFrac distance). For the taxonomy-level analysis, users can normalize the data using the widely used centered log-ratio (CLR) transformation method to relax the compositional constraint of the data, yet three other normalization methods of arcsine-root, rarefied count and proportion are also available. For reference, users can download all the resulting alpha and beta diversity indices and normalized taxonomic data.

**(4) Non-Mediational Analysis** - In this module, users can perform some baseline exploratory survival analysis to survey the disparity in survival time between medical treatments (e.g., treatment vs. placebo, new treatment vs. old treatment, and so forth) or environmental exposures (e.g., rural vs. urban, smoking vs. non-smoking, calorie restriction vs. ad libitum diet, and so forth). For this, users need to select (i) a survival time variable, (ii) a censored/event indicator variable, (iii) a treatment variable, (iv) covariate(s), and (v) an analytic method.

**(5) Mediational Analysis: Community-level Analysis: Alpha Diversity** - In this module, users can perform the microbiome causal mediation analysis to test jointly if the treatment alters the microbial alpha diversity, and then the altered microbial alpha diversity, in turn, influences the survival responses. For this, users need to select (i) a survival time variable, (ii) a censored/event indicator variable, (iii) a treatment variable, and (iv) covariate(s) for covariate-adjusted analysis or not for univariate analysis. The available analytic method here is the Imai method coupled with the Weibull regression model.

**(6) Mediational Analysis: Community-level Analysis: Beta Diversity** - In this module, users can perform the microbiome causal mediation analysis to test jointly if the treatment alters the microbial beta diversity, and then the altered microbial beta diversity, in turn, influences the survival responses. For this, users need to select (i) a survival time variable, (ii) a censored/event indicator variable, (iii) a treatment variable, and (iv) covariate(s) for covariate-adjusted analysis or not for univariate analysis. The available analytic method here is DACT coupled with MiRKAT for the treatment model and MiRKAT-S for the outcome model.

**(7) Mediational Analysis: Taxonomy-level Analysis** - In this module, users can perform the microbiome causal mediation analysis to test jointly if the treatment alters microbial taxa, and then the altered microbial taxa, in turn, influence the survival responses. For this, users need to select (i) a data format (default: CLR), (ii) a survival time variable, (iii) a censored/event indicator variable, (iv) a treatment variable, and (v) covariate(s) for covariate-adjusted analysis or not for univariate analysis, and (vi) the taxonomic ranks to be analyzed from phylum to genus for 16S ribosomal RNA gene amplicon sequencing or from phylum to species for shotgun metagenomics. The available analytic method here is the Imai method coupled with the Weibull regression model.

## GitHub Repository Contents

**(1) Data** - In this directory, example microbiome data are stored in a widely used unified format, called phyloseq, that can efficiently combine all essential microbiome data components as well as individual files.

**(2) Source** - In this directory, all the R functions that are needed to run MiMedSurv are stored.

**(3) www** - In this directory, some photos that are used to decorate the GUI of MiMedSurv are stored.

**(4) app.R** - In this file, all the central codes to control for user-interfaces and server functions of MiMedSurv are stored.

## Prerequites

#### Notice: For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. 

shiny
```
install.packages('shiny')
```

## Launch App

```
library(shiny)

runGitHub('MiMedSurvGit', 'yj7599', ref = 'main')
```

## Troubleshooting Tips

If you have any problems for using MiMedSurv, please report in issues (https://github.com/YJ7599/mimedsurvgit/issues) or email Hyo Jung Jang (hyojung.jang@northwestern.edu).





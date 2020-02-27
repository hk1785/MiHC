# R package: MiHC

Title: Microbiome Higher Criticism Analysis

Version: 1.0

Date: 2020-02-26

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hkoh7@jhu.edu>

Description: This R package (MiHC v1.0) provides facilities for MiHC which tests the association between a microbial group (e.g., community or clade) composition and a host phenotype of interest. MiHC is a data-driven omnibus test taken in a search space spanned by tailoring the higher criticism test to incorporate phylogenetic information and/or modulate sparsity levels and including the Simes test for excessively high sparsity levels. 

NeedsCompilation: No

Depends: R(>= 3.4.1)

Imports: cluster, compositions, permute, phyloseq

Suggests: knitr, rmarkdown

License: GPL-2

URL: https://github.com/hk1785/MiHC

## Reference

* Koh, H., Zhao, N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. (_under revision_).

## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/MiHC/issues) or email Hyunwook Koh (hkoh@jhu.edu).

* Tip 1. Depending on your pre-installed R libraries, this R package can require you to install additional R packages such as "gh", "usethis", "cli", etc using the command: install.packages("package_name").
* Tip 2. Please make sure if you have the most recent package version.

## Prerequites


cluster
```
install.packages("cluster")
```
compositions
```
install.packages("compositions")
```
devtools
```
install.packages("devtools")
```
permute
```
install.packages("permute")
```
phyloseq
```
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
```

## Installation

```
library(devtools)
install_github("hk1785/MiHC", force=T)
```

## Data format

```
library(phyloseq)
URL: https://joey711.github.io/phyloseq/
```

---------------------------------------------------------------------------------------------------------------------------------------

# Manual
This R package includes two core functions, MiHC and MiHC.plot. Please find the details below.

## :mag: MiHC

### Description
This function tests the association between a microbial group (e.g., community or clade) composition and a host phenotype of interest using MiHC.

### Usage
```
MiHC(y, covs, otu.tab, tree, model, hs=c(1,3,5,7,9), W=TRUE, comp=FALSE, CLR=FALSE, opt.ncl=30, n.perm=5000)
```

### Arguments
* _y_ - A numeric vector of the host outcomes. Gaussian (e.g., body mass index), Binomial (e.g., disease status, treatment/placebo) or Poisson (e.g., number of tumors/treatments) outcomes.
* _covs_ - A data.frame (or matrix/vector) for covariate (e.g., age, gender) adjustment(s). Default is cov=NULL for no covariate adjustment.
* _otu.tab_ - A matrix of the OTU table. (1. Rows are samples and columns are OTUs. 2. Monotone/singletone OTUs need to be removed.)
* _tree_ - A rooted phylogenetic tree.
* _model_ - "gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes, "poisson" for Poisson outcomes.
* _hs_ - A vector of the candidate modulation schema for lower sparsity levels. Default is hc=c(1,3,5,7,9).

* _W_ - An indicator to consider weighted high criticism tests or not. Default is W=TRUE to consider weighted higher criticism tests.
* _comp_ - An indicator if the OTU table contains absolute abundances (i.e., counts) or relative abundances (i.e., proportions). Default is comp=FALSE for absolute abundances.
* _CLR_ - An indicator if the OTU table needs to be converted using the centered log-ratio (CLR) transformation. Default is CLR=FALSE for no CLR transformation.
* _opt.ncl_ - A upper limit to find the optimal number of clusters. Default is opt.ncl=30.
* _n.perm_ - A number of permutations. Default is n.perm=5000. 

### Values
_$simes.pv_ - The p-value for the Simes test.

_$ind.pvs_ - The p-values for the item-by-item unweighted and weighted higher criticism tests.

_$ada.pvs_ - The p-values for the local (i.e., uHC(A) and wHC(A)) and global (i.e., MiHC) omnibus higher criticism tests.

### References
* Koh and Zhao. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. (_Under revision_).

* Simes (1986). An improved Bonferroni procedure for multiple tests of significance. _Biometrika_ 73(3):751-754

### Example
Import requisite R packages
```
library(cluster)
library(compositions)
library(permute)
library(phyloseq)
library(MiHC)
```
Import example microbiome data
```
data(phy)
otu.tab <- otu_table(phy)
tree <- phy_tree(phy)
y <- sample_data(phy)$y
covs <- data.frame(matrix(NA, length(y), 2))
covs[,1] <- as.numeric(sample_data(phy)$x1)
covs[,2] <- as.factor(sample_data(phy)$x2)
```

Fit MiHC
```
set.seed(123)
out <- MiHC(y, covs=covs, otu.tab=otu.tab, tree=tree, model="binomial")
out
```

## :mag: MiHC.plot

### Description
The Q-Q plots for the microbiome higher criticism analysis

### Usage
```
MiHC.plot(MiHC.out, leg.loc="bottomright", pdf.filename=NULL)
```

### Arguments
* _MiHC.out_ - An output obtained using the MiHC function.
* _leg.loc_ - The legend location to list the top 10 influential OTUs. Default is leg.loc="bottomright".
* _pdf.filename_ - The PDF filename to print the figure as a PDF file. Default is pdf.filename=NULL to print the figure on the R graphics window.

### Values
The Q-Q plots between the expected and observed quantiles for the unweighted and weighted higher criticism tests. Blue dots represent individual OTUs and a red diagonal line represents no influential points; as such, the OTUs that fall along the diagonal line have no influence on the host phenotype while the OTUs that have larger deviations from the diagonal line are more influential on the host phenotype. Darker to lighter vertical lines represent more to less influential OTUs in rank order among the 10 most influential OTUs that correspond to the 10 largest deviations from the red diagonal line.

### References
* Koh, H., Zhao, N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. (_under revision_).

### Example
Import requisite R packages
```
library(cluster)
library(compositions)
library(permute)
library(phyloseq)
library(MiHC)
```
Import example microbiome data
```
data(phy)
otu.tab <- otu_table(phy)
tree <- phy_tree(phy)
y <- sample_data(phy)$y
covs <- data.frame(matrix(NA, length(y), 2))
covs[,1] <- as.numeric(sample_data(phy)$x1)
covs[,2] <- as.factor(sample_data(phy)$x2)
```

Fit MiHC
```
set.seed(123)
out <- MiHC(y, covs=covs, otu.tab=otu.tab, tree=tree, model="binomial")
out
```

Create a graph
```
MiHC.plot(out)
```

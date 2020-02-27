\name{MiHC}
\alias{MiHC}
\title{
Microbiome higher criticism analysis
}
\description{
This function tests the association between a microbial group (e.g., community or clade) composition and a host phenotype of interest using MiHC.
}
\usage{
MiHC(y, covs=NULL, otu.tab, tree, model, hs=c(1,3,5,7,9), W=TRUE, 
comp=FALSE, CLR=FALSE, opt.ncl=30, n.perm=5000)
}
\arguments{
  \item{y}{
A numeric vector of the host outcomes. Gaussian (e.g., body mass index), Binomial (e.g., disease status, treatment/placebo) or Poisson (e.g., number of tumors/treatments) outcomes.
}
  \item{covs}{
A data.frame (or matrix/vector) for covariate (e.g., age, gender) adjustment(s). Default is cov=NULL for no covariate adjustment.
}
  \item{otu.tab}{
A matrix of the OTU table. (1. Rows are samples and columns are OTUs. 2. Monotone/singletone OTUs need to be removed.)
}
  \item{tree}{
A rooted phylogenetic tree.
}
  \item{model}{
"gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes, "poisson" for Poisson outcomes.
}
  \item{hs}{
A vector of the candidate modulation schema for lower sparsity levels. Default is hc=c(1,3,5,7,9).
}
  \item{W}{
An indicator to consider weighted high criticism tests or not. Default is W=TRUE to consider weighted higher criticism tests.
}
  \item{comp}{
An indicator if the OTU table contains absolute abundances (i.e., counts) or relative abundances (i.e., proportions). Default is comp=FALSE for absolute abundances.
}
  \item{CLR}{
An indicator if the OTU table needs to be converted using the centered log-ratio (CLR) transformation. Default is CLR=FALSE for no CLR transformation.
}
  \item{opt.ncl}{
A upper limit to find the optimal number of clusters. Default is opt.ncl=30.
}
  \item{n.perm}{
A number of permutations. Default is n.perm=5000.
}
}
\value{
simes.pv: The p-value for the Simes test.

ind.pvs: The p-values for the item-by-item unweighted and weighted higher criticism tests.

ada.pvs: The p-values for the local (i.e., uHC(A) and wHC(A)) and global (i.e., MiHC) omnibus higher criticism tests.

}
\references{
Koh and Zhao. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. (Under revision).

Simes (1986). An improved Bonferroni procedure for multiple tests of significance. Biometrika.73(3):751-754

}
\author{
Hyunwook Koh
}

\examples{
# Import requisite R packages
require(cluster)
require(compositions)
require(permute)
require(phyloseq)

# Import example microbiome data
data(phy)
otu.tab <- otu_table(phy)
tree <- phy_tree(phy)
y <- sample_data(phy)$y
covs <- data.frame(matrix(NA, length(y), 2))
covs[,1] <- as.numeric(sample_data(phy)$x1)
covs[,2] <- as.factor(sample_data(phy)$x2)

# Fit MiHC
set.seed(123)
out <- MiHC(y, covs=covs, otu.tab=otu.tab, tree=tree, model="binomial", n.perm=1000)
out

}

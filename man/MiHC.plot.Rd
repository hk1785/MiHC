\name{MiHC.plot}
\alias{MiHC.plot}
\title{
The Q-Q plots for the microbiome higher criticism analysis
}
\description{
This function produces Q-Q plots for the microbiome higher criticism analysis.
}
\usage{
MiHC.plot(MiHC.out, leg.loc="bottomright", pdf.filename=NULL)
}
\arguments{
  \item{MiHC.out}{
An output obtained using the MiHC function.
}
  \item{leg.loc}{
The legend location to list the top 10 influential OTUs. Default is leg.loc="bottomright".
}
  \item{pdf.filename}{
The PDF filename to print the figure as a PDF file. Default is pdf.filename=NULL to print the figure on the R graphics window.
}
}
\value{
The Q-Q plots between the expected and observed quantiles for the unweighted and weighted higher criticism tests. Blue dots represent individual OTUs and a red diagonal line represents no influential points; as such, the OTUs that fall along the diagonal line have no influence on the host phenotype while the OTUs that have larger deviations from the diagonal line are more influential on the host phenotype. Darker to lighter vertical lines represent more to less influential OTUs in rank order among the 10 most influential OTUs that correspond to the 10 largest deviations from the red diagonal line.

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

# Create plots
MiHC.plot(MiHC.out=out)

}

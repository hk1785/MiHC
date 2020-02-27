\name{MiHC.stat}
\alias{MiHC.stat}
\title{
uHC(h)'s and wHC(h)'s
}
\description{
This function produces the values of the test statistics of uHC(h)'s and wHC(h)'s.
}
\usage{
MiHC.stat(Zs, hs=hs, Ws=Ws)
}
\arguments{
  \item{Zs}{
The Z-score values for individual OTUs.
}
  \item{hs}{
A vector of the candidate modulation schema for lower sparsity levels. Default is hc=c(1,3,5,7,9).
}
  \item{Ws}{
A vector of the weights for invidual OTUs.
}
}
\value{
The values of the test statistics of uHC(h)'s and wHC(h)'s.

}
\references{
Koh and Zhao. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. (Under revision).

}
\author{
Hyunwook Koh
}

\examples{
# Import requisite R packages
require(cluster)
require(permute)
require(phyloseq)

Zs <- rnorm(50)
hs <- c(1,3,5,7,9)
Ws <- runif(50)

MiHC.stat(Zs, hs=hs, Ws=Ws)

}
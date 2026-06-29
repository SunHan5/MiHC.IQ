\name{HC.IQ}
\alias{HC.IQ}
\title{Higher Criticism Test with Integrated Quantiles}

\description{
  This function implements the Higher Criticism test with integrated quantile scores for microbiome association analysis. It constructs sample weights based on quantile regression (either integrated over all quantiles or at a specific quantile), computes marginal association statistics for each OTU, and performs modulated higher criticism tests with optional phylogenetic weighting. The method is designed to detect sparse association signals while accommodating non-normal and heterogeneous phenotype distributions.
}

\usage{
HC.IQ(
  y,
  covs = NULL,
  otu.tab,
  tree,
  tau = -1,
  score,
  hs = c(1, 3, 5, 7, 9),
  W = TRUE,
  comp = FALSE,
  CLR = FALSE,
  opt.ncl = 30,
  n.perm = 5000,
  seed = 123
)
}

\arguments{
  \item{y}{
    A numeric vector of length \eqn{n} representing the host phenotype (continuous response variable).
  }
  \item{covs}{
    An optional data.frame, matrix, or vector of covariates (e.g., age, sex). Default is \code{NULL}.
  }
  \item{otu.tab}{
    An \eqn{n \times p} matrix or data.frame of OTU abundances. Rows correspond to samples and columns to OTUs. The function optionally performs composition-normalization and CLR transformation.
  }
  \item{tree}{
    A phylogenetic tree of class \code{phylo} (from the \code{ape} package) describing the evolutionary relationships among OTUs. Used for computing cophenetic distances and phylogenetic weights.
  }
  \item{tau}{
    Quantile level. If \code{tau = -1} (default), the integrated quantile score is constructed using rank integration across all quantiles. If \code{tau} is set to a specific value in \eqn{(0,1)}, the score is constructed using the sub-gradient of quantile regression at that quantile.
  }
  \item{score}{
    A character string specifying the quantile weight function. Options are: \code{"wilcoxon"}, \code{"normal"}, \code{"lehmann"}, or \code{"inverselehmann"}. This argument is used when \code{tau = -1}.
  }
  \item{hs}{
    A numeric vector of truncation parameters \eqn{h} for the modulated higher criticism statistic. Default is \code{c(1, 3, 5, 7, 9)}.
  }
  \item{W}{
    A logical value. If \code{TRUE} (default), phylogenetic weighting is applied to the higher criticism statistic. If \code{FALSE}, only unweighted version is computed.
  }
  \item{comp}{
    A logical value. If \code{FALSE} (default), OTU abundances are normalized to relative abundances (each row sums to 1). If \code{TRUE}, the input \code{otu.tab} is assumed to be already compositional.
  }
  \item{CLR}{
    A logical value. If \code{TRUE}, the centered log-ratio (CLR) transformation is applied to the OTU abundance matrix. Default is \code{FALSE}.
  }
  \item{opt.ncl}{
    An integer specifying the maximum number of clusters for PAM (Partitioning Around Medoids) clustering of OTUs based on phylogenetic distances. The optimal number is determined by the maximum average silhouette width. Default is \code{30}. Special values \code{"gmx"} and \code{"fmx"} are also supported for automatic cluster selection using specific algorithms.
  }
  \item{n.perm}{
    Number of permutations for constructing the null distribution. Default is \code{5000}.
  }
  \item{seed}{
    A numeric seed for reproducible permutation results. Default is \code{123}.
  }
}

\value{
  A list containing the following components:
  \describe{
    \item{graph.els}{A list of elements for visualization:
      \describe{
        \item{otu.ids}{Character vector of OTU identifiers.}
        \item{exp.HC}{Expected higher criticism values under the null.}
        \item{obs.HC}{Observed higher criticism values.}
        \item{exp.wHC}{Expected weighted higher criticism values (if \code{W = TRUE}).}
        \item{obs.wHC}{Observed weighted higher criticism values (if \code{W = TRUE}).}
      }
    }
    \item{opt.nclust}{Optimal number of clusters selected by PAM (if \code{W = TRUE}).}
    \item{ind.pvs}{A matrix of p-values for individual HC tests:
      \itemize{
        \item For \code{W = TRUE}: columns \code{"uHC(h)"} and \code{"wHC(h)"} for each truncation parameter \eqn{h}.
        \item For \code{W = FALSE}: a single column \code{"uHC(h)"}.
      }
    }
    \item{ada.pvs}{A vector of adaptive p-values:
      \itemize{
        \item For \code{W = TRUE}: \code{"uHC(A)"} (unweighted adaptive), \code{"wHC(A)"} (weighted adaptive), and \code{"MiHC"} (omnibus).
        \item For \code{W = FALSE}: \code{"uHC(A)"} only.
      }
    }
  }
}

\details{
  The function proceeds as follows:

  \enumerate{
    \item \strong{Sample weight construction:}
      If \code{tau = -1}, the integrated quantile score is computed using the rank integration method (via \code{MiRKAT:::NewRanks}). If \code{tau} is a specific quantile, weights are constructed using the sub-gradient of quantile regression residuals.

    \item \strong{Marginal association statistics:}
      For each OTU, the raw association score \eqn{S_j = \sum_i w_i X_{ij}} is computed. Permutation is used to estimate the standard error, yielding standardized Z-statistics \eqn{Z_j}.

    \item \strong{Modulated higher criticism:}
      Marginal p-values are computed from \eqn{Z_j}. The higher criticism statistic is computed for each truncation parameter \eqn{h \in} \code{hs}.

    \item \strong{Phylogenetic weighting (if \code{W = TRUE}):}
      Cophenetic distances are computed from the phylogenetic tree. OTUs are clustered using PAM, and phylogenetic weights are computed based on within-cluster distances and marginal effect sizes.

    \item \strong{Adaptive integration:}
      A double-permutation strategy is used to combine individual tests and compute omnibus p-values.
  }
}

\examples{
# Load required packages
library(phyloseq)

# Load CRC dataset
data(CRC)

# Prepare data
otu.tab <- otu_table(CRC)
tree <- phy_tree(CRC)
sample.data <- sample_data(CRC)
y <- as.numeric(sample.data$BMI)
cov <- data.frame(Age = as.numeric(sample.data$Age))

# Run HC.IQ with integrated quantile score (Wilcoxon)
result <- HC.IQ(
  y = y,
  covs = cov,
  otu.tab = otu.tab,
  tree = tree,
  tau = -1,
  score = "wilcoxon",
  hs = c(1, 3, 5, 7, 9),
  W = TRUE,
  n.perm = 5000
)

# Examine results
result$ada.pvs
result$ind.pvs
}

\references{
  Koh H, Zhao N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. \emph{Microbiome}. 2020;8(1):63.

  Wang T, Ling W, Plantinga AM, Wu MC, Zhan X. Testing microbiome association using integrated quantile regression models. \emph{Bioinformatics}. 2022;38(2):419-425.
}

\seealso{
  \code{\link{MiHC.IQ}} for the omnibus wrapper function that combines results across multiple quantile weight functions.
}
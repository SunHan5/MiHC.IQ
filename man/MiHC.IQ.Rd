\name{MiHC.IQ}
\alias{MiHC.IQ}
\title{Omnibus Higher Criticism Test with Integrated Quantiles}

\description{
  This function serves as an omnibus wrapper for \code{\link{HC.IQ}} that simultaneously evaluates multiple quantile weighting functions. It applies the Higher Criticism test with integrated quantile scores using each specified weighting function and combines the results using ACAT (Aggregated Cauchy Association Test) to produce an omnibus test. This approach automatically selects the optimal weighting strategy, enhancing robustness across different distributional patterns of the host phenotype.
}

\usage{
MiHC.IQ(
  y,
  covs,
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
    A data.frame, matrix, or vector of covariates (e.g., age, sex).
  }
  \item{otu.tab}{
    An \eqn{n \times p} matrix or data.frame of OTU abundances. Rows correspond to samples and columns to OTUs.
  }
  \item{tree}{
    A phylogenetic tree of class \code{phylo} (from the \code{ape} package) describing the evolutionary relationships among OTUs.
  }
  \item{tau}{
    Quantile level. If \code{tau = -1} (default), the integrated quantile score is constructed using rank integration across all quantiles. If \code{tau} is set to a specific value in \eqn{(0,1)}, the score is constructed using the sub-gradient of quantile regression at that quantile.
  }
  \item{score}{
    A character vector specifying one or more quantile weighting functions to be evaluated. Options are: \code{"wilcoxon"}, \code{"normal"}, \code{"lehmann"}, and \code{"inverselehmann"}. The omnibus test automatically combines results across all specified weights using ACAT.
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
    An integer specifying the maximum number of clusters for PAM (Partitioning Around Medoids) clustering of OTUs based on phylogenetic distances. The optimal number is determined by the maximum average silhouette width. Default is \code{30}.
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
    \item{ind.pvs}{A matrix of p-values for individual HC tests. Rows are named as \code{"uHC1"}, \code{"uHC3"}, ..., \code{"wHC9"} (unweighted and weighted HC for each truncation parameter \eqn{h}). Columns correspond to the weighting functions specified in \code{score} (e.g., \code{"Wilcoxon"}, \code{"Normal"}, \code{"Lehmann"}, \code{"Inverselehmann"}).}
    \item{ada.pvs}{A named vector of adaptive p-values:
      \itemize{
        \item \code{"MiHC-IQW"}: omnibus p-value for the Wilcoxon weighting.
        \item \code{"MiHC-IQN"}: omnibus p-value for the Normal weighting.
        \item \code{"MiHC-IQI"}: omnibus p-value for the Inverse-Lehmann weighting.
        \item \code{"MiHC-IQL"}: omnibus p-value for the Lehmann weighting.
        \item \code{"MiHC-IQ"}: optimal p-value combining all weighting functions via ACAT.
      }
      Only the weighting functions specified in \code{score} are returned.
    }
  }
}

\details{
  The \code{MiHC.IQ} function iterates over the specified weighting functions in \code{score}, calls \code{\link{HC.IQ}} for each, and collects the results. The omnibus test for each weighting function selects the minimum p-value across all truncation parameters and phylogenetic weighting options (uHC vs. wHC).

  To combine evidence across different weighting functions, the function applies the Aggregated Cauchy Association Test (ACAT) to the omnibus p-values from each weighting function. This produces the final \code{"MiHC-IQ"} p-value, which is robust to varying quantile shift patterns without requiring prior specification of the optimal weighting function.

  For details on the underlying HC-IQ methodology, please refer to \code{\link{HC.IQ}}.
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

# Run MiHC.IQ with all four weighting functions
out <- MiHC.IQ(
  y = y,
  covs = cov,
  otu.tab = otu.tab,
  tree = tree,
  tau = -1,
  score = c("wilcoxon", "normal", "lehmann", "inverselehmann"),
  hs = c(1, 3, 5, 7, 9),
  W = TRUE,
  n.perm = 5000
)

# Extract the individual p-value matrix
out$ind.pvs

# Examine the adaptive p-values
out$ada.pvs
}

\references{
  Wang T, Ling W, Plantinga AM, Wu MC, Zhan X. Testing microbiome association using integrated quantile regression models. \emph{Bioinformatics}. 2022;38(2):419-425.

  Liu Y, Chen S, Li Z, et al. ACAT: a fast and powerful p value combination method for rare-variant analysis in sequencing studies. \emph{The American Journal of Human Genetics}. 2019;104(3):410-421.

  Koh H, Zhao N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. \emph{Microbiome}. 2020;8(1):63.
}

\seealso{
  \code{\link{HC.IQ}} for the underlying function that performs the test for a single weighting function.
}
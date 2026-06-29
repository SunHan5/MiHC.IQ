# MiHC.IQ

**Type:** Package

**Title:** Microbiome higher criticism analysis with integrated quantile regression models

**Version:** 1.0

**Author:** Han Sun, Peilin Pang, Yuejiao Qiao, Jinglei Song, Jianjun Zhang, Jianping Wang, Keqin Su, Liangliang Liu, Yuying Chen, Yuting Tan, Ban Huo, Ge Song

**Maintainer:** Han Sun <sunh529@henau.edu.cn>

**Imports:** phyloseq, compositions, quantreg, permute, cluster, MiRKAT, MiHC, ACAT

**Description:** MiHC.IQ implements a unified analytical framework that integrates quantile regression with modulated higher criticism statistics for microbiome association analysis. It simultaneously addresses signal sparsity and distributional heterogeneity of host phenotypes, offering robust detection of sparse microbial association signals across the entire phenotypic distribution.

**License:** GPL-2

**Encoding:** UTF-8

**LazyData:** true

**URL:** https://github.com/SunHan5/MiHC.IQ


## Introduction

**MiHC.IQ** is an R package that implements a unified analytical framework integrating quantile regression with modulated higher criticism statistics for microbiome association analysis. It simultaneously addresses:

- **Detect sparse signals**: using higher criticism statistics to aggregate marginal signals across OTUs, thereby performing a community-level association test while remaining sensitive to sparse signals.
- **Capture distributional heterogeneity**: using integrated quantile scores that encode the entire phenotypic distribution into sample weights.
- **Leverage phylogenetic information**: using phylogenetically weighted HC statistics to incorporate evolutionary relationships among OTUs.

The method applies a double-permutation strategy to rigorously control type I error while adaptively combining candidate tests across different truncation parameters and weighting schemes (Wilcoxon, Normal, Lehmann, and Inverse-Lehmann). The final optimal p-value is obtained via ACAT (Aggregated Cauchy Association Test).

**MiHC.IQ** can be applied to cross-sectional microbiome data with continuous host phenotypes (e.g., BMI, blood pressure).


## Installation

You may install `compositions, quantreg, permute, cluster, MiRKAT` packages through the following code:

```
install.packages(c("compositions", "quantreg", "permute", "cluster", "MiRKAT"))
```

```
BiocManager::install("phyloseq")
```

```
devtools::install_github("yaowuliu/ACAT")
```

You may install `MiHC.IQ` from GitHub using the following code: 

```
devtools::install_github("SunHan5/MiHC.IQ", force = TRUE)
```


## Usage
```
MiHC.IQ(y, covs, otu.tab, tree, tau = -1, score = c("wilcoxon", "normal", "lehmann", "inverselehmann"), hs = c(1, 3, 5, 7, 9), W = TRUE,  comp = FALSE, CLR = FALSE, opt.ncl = 30, n.perm = 5000, seed = 123)
```


## Arguments
* _y_ - Numeric vector of host phenotype (continuous).
* _covs_ - covariate (e.g., age, gender). Default is covs = NULL.
* _otu.tab_ - A matrix of the OTU table. (1. Rows are samples and columns are OTUs. 2. Monotone/singletone OTUs need to be removed.)
* _tree_ - A rooted phylogenetic tree. 
* _tau_ - Quantile level (-1 for integrated, or value in (0,1)).
* _score_ - Weight function(s): `"wilcoxon"`, `"normal"`, `"lehmann"`, `"inverselehmann"`.
* _hs_ - Modulation parameters, default `c(1, 3, 5, 7, 9)`.
* _W_ - Logical, whether to apply phylogenetic weighting.
* _comp_ - An indicator if the OTU table contains absolute abundances or relative abundances. Default is comp = FALSE for absolute abundances.
* _CLR_ - An indicator if the OTU table needs to be converted using the centered log-ratio (CLR) transformation. Default is CLR=FALSE for no CLR transformation.
* _opt.ncl_ - A upper limit to find the optimal number of clusters. Default is opt.ncl=30.
* _n.perm_ - A number of permutations. Default is n.perm = 5000. 
* _seed_ - Random number generator. Default is seed = 123.


## Values
_$Ind.pvs_: The p-values for individual tests.

_$ada.pvs_: The p-values for the omnibus tests and optimal test.


## Example

Import requisite R packages: 

```
library(phyloseq)
library(compositions)
library(quantreg)
library(permute)
library(cluster)
library(MiRKAT)
library(ACAT)
library(MiHC.IQ)
```


Import example microbiome data:

```
data(CRC)
otu.tab <- otu_table(CRC)
tree <- phy_tree(CRC)
sample.data <- sample_data(CRC)
y <- as.numeric(sample.data$BMI)
cov <- data.frame(Age = as.numeric(sample.data$Age))
```


Fit MiHC.IQ:

```
set.seed(123)
out <- MiHC.IQ(y = y, covs = cov, otu.tab = otu.tab, tree = tree, score = c("wilcoxon", "normal", "inverselehmann", "lehmann"))
out
```


## References

* Koh H, Zhao N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. Microbiome. 2020;8(1):63.

* Wang T, Ling W, Plantinga AM, et al. Testing microbiome association using integrated quantile regression models. Bioinformatics. 2022;38(2):419-425.

* Liu Y, Chen S, Li Z, et al. ACAT: a fast and powerful p value combination method for rare-variant analysis in sequencing studies. The American Journal of Human Genetics. 2019;104(3):410-421.

* McMurdie PJ, Holmes S. phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE. 2013;8(4):e61217.

* Baxter NT, Ruffin MT, Rogers MAM, et al. Microbiota-based model improves the sensitivity of fecal immunochemical test for detecting colonic lesions. Genome Medicine. 2016;8(1):37.


## Statement

Our code mainly refers to R packages, _MiHC_ and _MiRKAT_.# MiHC.IQ

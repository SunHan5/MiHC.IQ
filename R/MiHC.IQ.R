#' @import compositions quantreg permute cluster MiRKAT MiHC ACAT
#' @export
MiHC.IQ <- function (y, covs, otu.tab, tree, tau = -1, score, hs = c(1, 3, 5, 7, 9), W = TRUE, comp = FALSE, CLR = FALSE, opt.ncl = 30, n.perm = 5000, seed = 123) {
  ind.pvs <- ada.pvs <- c()
  
  if ("wilcoxon" %in% score) {
    out.wilcoxon <- HC.IQ(y, covs, otu.tab, tree, tau, score = "wilcoxon")
    ind.pvs <- cbind(ind.pvs, c(out.wilcoxon$ind.pvs))
    ada.pvs <- c(ada.pvs, out.wilcoxon$ada.pvs["MiHC"])
    colnames(ind.pvs)[ncol(ind.pvs)] <- "Wilcoxon"
    names(ada.pvs)[length(ada.pvs)] <- "MiHC-IQW"
  }
  
  if ("normal" %in% score) {
    out.normal <- HC.IQ(y, covs, otu.tab, tree, tau, score = "normal")
    ind.pvs <- cbind(ind.pvs, c(out.normal$ind.pvs))
    ada.pvs <- c(ada.pvs, out.normal$ada.pvs["MiHC"])
    colnames(ind.pvs)[ncol(ind.pvs)] <- "Normal"
    names(ada.pvs)[length(ada.pvs)] <- "MiHC-IQN"
  }
  
  if ("inverselehmann" %in% score) {
    out.inverselehmann <- HC.IQ(y, covs, otu.tab, tree, tau, score = "inverselehmann")
    ind.pvs <- cbind(ind.pvs, c(out.inverselehmann$ind.pvs))
    ada.pvs <- c(ada.pvs, out.inverselehmann$ada.pvs["MiHC"])
    colnames(ind.pvs)[ncol(ind.pvs)] <- "Inverselehmann"
    names(ada.pvs)[length(ada.pvs)] <- "MiHC-IQI"
  }
  
  if ("lehmann" %in% score) {
    out.lehmann <- HC.IQ(y, covs, otu.tab, tree, tau, score = "lehmann")
    ind.pvs <- cbind(ind.pvs, c(out.lehmann$ind.pvs))
    ada.pvs <- c(ada.pvs, out.lehmann$ada.pvs["MiHC"])
    colnames(ind.pvs)[ncol(ind.pvs)] <- "Lehmann"
    names(ada.pvs)[length(ada.pvs)] <- "MiHC-IQL"
  }

  ada.pvs <- c(ada.pvs, ACAT(ada.pvs))
  rownames(ind.pvs) <- c(paste("uHC", seq(1, 9, 2), sep = ""), paste("wHC", seq(1, 9, 2), sep = ""))
  names(ada.pvs)[length(ada.pvs)] <- "MiHC-IQ"
  return(list(ind.pvs = ind.pvs, ada.pvs = ada.pvs))
}

#' Calculates the logloss of model predictions used in caret package
#'
#' @param data data as in caret defined
#' @param lev levels of the classes used for classifying
#' @param model not used
#'
#' @return
#' @export
#'
#' @examples
LogLoss = function(data, lev=NULL, model=NULL) { 
  epp = 10^-15
  padj = pmax(pmin(data[,4],1-epp),epp)
  yi = rep(0,nrow(data))
  yi[which(data$obs == lev[2])] = 1
  ll = -mean(yi * log(padj) + (1 - yi)* log(1-padj), na.rm = T)
  names(ll) = c("logLoss")
  ll
}
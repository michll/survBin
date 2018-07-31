#' Generates stratified folds used for cross-validation
#'
#' @param y binary variable eg. (event, no event)
#' @param k number of folds (default: 10)
#' @param stratified stratification (default: TRUE)
#'
#' @return indices of cross-validation folds
#' @export
#'
#' @examples
crossvalFolds <- function(y, k=10, stratified=T) {
  if (stratified) {
    t <- table(y)
    f <- matrix(nrow=length(t), ncol=k)
    for (i in 1:length(t)) {
      f[i,] <- ipred::kfoldcv(k, t[i])
    }
    rownames(f) <- names(t)
    idx <- rep(NA, length(y))
    data <- data.frame(y=y, x=1:length(y))
    for (i in 1:(k-1)) {
      s <- sampling::strata(data, "y", f[match(unique(data$y), rownames(f)), i], "srswor")
      sy <- table(data[s$ID, "y"])
      stopifnot(sy[rownames(f)] == f[,i])
      idx[data[s$ID, "x"]] <- i
      stopifnot(sum(!is.na(idx)) == sum(f[, 1:i]))
      data <- data[-s$ID,]
    }
    idx[data[, "x"]] <- k
    stopifnot(all(!is.na(idx)))
    return (idx)
  } else {
    folds.size <- kfoldcv(k, length(y))
    idx <- c()
    for (i in 1:k) {
      idx <- c(idx, rep(i, folds.size[i]))
    }
    idx <- sample(idx)
    return (idx)
  }
}

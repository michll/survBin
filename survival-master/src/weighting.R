weighting = function(x, sq = 0.9){
  xnew = x * 3/4 * (1-seq(-sq, sq, length.out = length(x))^2)
  return(xnew)
}
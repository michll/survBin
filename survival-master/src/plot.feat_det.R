plot.feat_det = function(x, scale_x = 1, colused = NULL, ncut=75, spty = 0.5){
  
  if(is.null(colused)){col1 = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6","#798E87", "#C27D38", "#CCC591", "#29211F")}
  else{col1 = colused}
  #parold = par()
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  ii = which(apply(x$vi_mat,1,max, na.rm = T)>ncut)
  imax = which(names(x$imp)[1] == rownames(x$vi_mat))
  par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE )
  plot(x$tseq/scale_x,x$vi_mat[imax,], type = "l", lwd=0.15,lty = 3, col = add.alpha("chartreuse3",1), xlab = "Time from start", 
       ylab = "Feature importance", ylim = c(0,100))
  lines(smooth.spline(x$tseq/scale_x,x$vi_mat[imax,], spar = spty), lwd=2, col = "chartreuse3")
  for(i in 2:min(length(x$ind),10)){
    inext = which(names(x$imp)[i] == rownames(x$vi_mat))
    lines(x$tseq/scale_x,x$vi_mat[inext,], type = "l", lwd=0.15, col = add.alpha(col1[i-1],1), lty=3)
    lines(smooth.spline(x$tseq/scale_x,x$vi_mat[inext,], spar = spty), lwd=2, col = col1[i-1])
  }
  legend(x = max(x$tseq/scale_x)*1.05,y = 110, bty = "n",legend = substr(names(x$imp)[1:min(length(x$ind),10)], start=1, stop=7), 
         col = c("chartreuse3", col1), lwd=2, lty=1, inset=c(-0.2,0))
  
}


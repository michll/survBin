#' Simulates a survival data set 
#'
#' @param nsub number of subjects
#' @param beta_vec vector for associations (main effects) between predictor and response
#' @param nfeat number of features (covariates)
#' @param correl correlation between features
#' @param lambda_vec used for generating the censoring distribution (baseline hazard)
#' @param surv_typ "cox" no other option used
#' @param ia effect of an interaction between first two main effects
#'
#' @return survival data set as a list with $x the predictor variables and $y a survival object (form Surv()) as the response
#' @export
#'
#' @examples
sim_surv = function(nsub,beta_vec,nfeat,correl,lambda_vec,surv_typ="cox", ia = NA){
  # simuliert cox Zufallsgroessen
  res=list()
  sigma=matrix(correl,nrow=nfeat,ncol=nfeat);diag(sigma)=1
  rr=mvtnorm::rmvnorm(nsub,mean=rep(0,nfeat),sigma=sigma)
  colnames(rr)=c(paste("x",1:length(beta_vec),sep=""),paste("r",1:(nfeat-length(beta_vec)),sep = ""))
  U1=runif(nsub,0,1);U2=runif(nsub,0,1)
  XB = rr[,1:length(beta_vec)]%*%beta_vec
  if(!is.na(ia)){
    XB = XB + rr[,1]*rr[,2]*ia 
  }
  Ti = -log(U1)/(lambda_vec[1]*exp(XB))
  Tc = -log(U2)/lambda_vec[2]
  To = pmin(Ti,Tc)  #observed time is min of censored and true
  event = To==Ti   # set to 1 if event is observed
  res$x = data.frame(rr)
  res$y = Surv(time = To+1,event = event)
  class(res$y) = "Surv"
  return(res)
}

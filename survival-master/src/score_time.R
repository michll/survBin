#' Calculates the time dependent AUC for a given risk score using the Simpsons rule for integration
#'
#' @param time event times, event time or censoring time
#' @param death (death) events, binary vector (0 for censored, and 1 for observed events)
#' @param riskScore risk score under evaluation
#' @param time_int time interval for prediction of time dependent AUC
#'
#' @return integrated AUC
#' @export
#'
#' @examples
score_time<-function(time, death, riskScore, time_int = c(0.1,5,0.1) ){
  # compute iAUC 
  times <- seq(time_int[1],time_int[2],by=time_int[3]) 
  aucs <- timeROC::timeROC(T=time,
                           delta=death,
                           marker=riskScore, 
                           cause=1,
                           weighting="marginal",
                           times=times,
                           iid=FALSE)$AUC
  
  # Simpsons rules for integrating under curve
  iAUC <- Bolstad2::sintegral(times, aucs)$int / (max(times) - min(times))
  return (iAUC)
}

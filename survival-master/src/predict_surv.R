#' Generates time dependent predictions from a surv_fct output
#'
#' @param surv_fct_model output of surv_fct
#' @param newx new data for prediction
#'
#' @return time dependent predictions 
#' @export
#'
#' @examples
predict_surv = function(surv_fct_model, newx){
  tsq = surv_fct_model$tseq
  datx = newx
  res = list()
  pred = matrix(NA, nrow = nrow(datx), ncol = length(tsq))
  for(i in 1:length(tsq)){
    
    m1 = surv_fct_model$model[[i]]
    if(is.na(m1)){
      pred[,i] = rep(NA, nrow(newx))
    }
    else{
      pred[,i] = predict(m1, newdata = as.matrix(datx), type = "prob")[,2]
    }
    # datx = cbind(datx, pred[,1:i])
  }
  res$pred = pred
  res$pred_mean = apply(pred,1,mean, na.rm = T)
  res$pred_int = apply(pred, 1, function(x) {
    ii = !is.na(x)
    Bolstad2::sintegral(tsq[ii], x[ii])$int / (max(tsq[ii]) - min(tsq[ii]))})
  return(res)
}

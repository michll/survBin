#' Generates a time dependent cox model matrix
#'
#' @param tsq 
#' @param m 
#' @param y 
#'
#' @return 
#' @export
#'
#' @examples
time_rs_fct = function(tsq, m, y){
  data_new = NULL
  for(i in 1:nrow(m)){
    cnew = NULL
    for(j in 1:ncol(m)){
      # id time event rs_t tstart tstop
      if(j == 1){
        cnew = rbind(cnew,c(paste0("id",i), y[i,1],y[i,2],ifelse(tsq[j]<y[i,1],0,ifelse(y[i,2]==1,1,NA)),  m[i,j], 0, tsq[j]))
      }
      else{
        cnew = rbind(cnew,c(paste0("id",i), y[i,1],y[i,2],ifelse(tsq[j]<y[i,1],0,ifelse(y[i,2]==1,1,NA)), m[i,j], tsq[j-1], tsq[j]))
      }
    }
    data_new = rbind(data_new, cnew)
  }
  data_new = data.frame(id = data_new[,1], time = as.numeric(data_new[,2]), event = as.numeric(data_new[,3]),
                        event_new = as.numeric(data_new[,4]),
                        rs_t = as.numeric(data_new[,5]), tstart = as.numeric(data_new[,6]), tstop = as.numeric(data_new[,7]))
  return(data_new)
}

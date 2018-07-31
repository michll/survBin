#' Basic function for evaluating different data sets or simulation setups
#'
#' @param data_input data input for prediction
#' @param k.cv number of cross-validation folds (default: 10)
#' @param seeed should the seed be fixed (default: 123)
#' @param data_src data source, one of "dream" for dream data input, "sim" for data simulation, "tcga" for TCGA data, or "other"
#' @param data_spec specification for data simulation in "sim" setup
#' @param ncl number of cores used for parallelization
#' @param tint time sequence for estimating the time dependent AUC curve, defined as c(start, stop, increments)
#'
#' @return prediction accuracy (time dependent AUC) per method
#' @export
#'
#' @examples sim_output = data_sel_eval(data_input, k.cv = 10, seeed = 123, data_src = "sim", 
#'                                      data_spec = list(n = c(200),
#'                                      p = c(50),
#'                                      corel = c(0),
#'                                      beta_vec = c(1,0.5,-1,-0.5),
#'                                      beta_ia = 1,
#'                                      surv_cens = c(0.5,0.55),
#'                                      rept = 2), 
#'                                      ncl = 2, tint = NULL)
data_sel_eval = function(data_input, k.cv = 10, seeed = 123, data_src = "dream", data_spec = NULL, ncl = 1, tint = NULL, dat1 = NULL){
  
  if(data_src=="dream"){
    not_to_use <- c("RPT", "STUDYID", "LKADT_P", "DEATH", "DISCONT", "ENDTRS_C", "ENTRT_PC")
    not_to_use_Idx <- which(colnames(data_input) %in% not_to_use)
    
    x <- model.matrix( ~ 0 + ., (data_input[, -not_to_use_Idx])[, apply(data_input[, -not_to_use_Idx], 2, function(l) !any(is.na(l)))])
    
    dat = res = list()
    dat$y = survival::Surv(time = data_input$LKADT_P, event = data_input$DEATH)
    dat$x = x
    dat1 = dat
    dat1$x = data_input[, -not_to_use_Idx]
    # remove columns with complete NAs
    dat1$x = dat1$x[,apply(dat1$x, 2, function(x) !anyNA(x))]
    
    # Function call to sim
    out = eval_data_cv(dat = dat, dat1 = dat1, k.cv = k.cv, seeed = seeed, time_inval = c(183,915,30.5), ncl = ncl)
  }
  if(data_src == "sim"){
    if(is.null(data_spec)){
      return("Provide data specifications for simulation setup.")
    }
    else{
      n = data_spec$n
      p = data_spec$p
      corel = data_spec$corel
      rept = data_spec$rept
      
      ll = expand.grid(n, p, corel)
      ll_rep = NULL
      for(i in 1:nrow(ll)){
        ll_rep = rbind(ll_rep, matrix(rep(ll[i,],rept), ncol = 3, byrow = T))
      }
      dat = list()
      set.seed(seeed)
      for(i in 1:nrow(ll_rep)){
        pj = as.numeric(ll_rep[i,2]); nj = as.numeric(ll_rep[i,1]); corj = as.numeric(ll_rep[i,3])
        dat[[i]] = sim_surv(nj, data_spec$beta_vec, pj, correl = corj, data_spec$surv_cens, ia = data_spec$beta_ia)
      }
    }
    
    out = lapply(dat, function(x) eval_data_cv(dat = x, 
                                               dat1 = x, 
                                               k.cv = k.cv, 
                                               seeed = seeed, 
                                               time_inval = c(round(quantile(x$y[,1],0.1),1),round(quantile(x$y[,1],0.9),1),0.1),
                                               ncl = ncl))
  }
  if(data_src == "tcga"){
    age = pData(TCGA_mutations)[,3]
    age[is.na(age)] = median(age, na.rm = T)
    mm = cbind(model.matrix( ~ ., pData(TCGA_mutations)[,c(4:5)])[,-1],age)
    ii = which(pData(TCGA_mutations)[,1]>0)
    ij = apply(t(exprs(TCGA_mutations)),2, function(x) ifelse(min(table(x))/length(x) < 0.015, FALSE, TRUE ))
    dat = list()
    dat$x = cbind(mm,t(exprs(TCGA_mutations))[,ij])[ii,]
    dat$y = survival::Surv(time = pData(TCGA_mutations)[,1], event = pData(TCGA_mutations)[,2])[ii,]
    dat1 = dat
    dat1$x = data.frame(age = age, gender = as.factor(pData(TCGA_mutations)[,4]), tumer_type = as.factor(pData(TCGA_mutations)[,5]),
                        t(exprs(TCGA_mutations))[,ij])[ii,]
    
    out = eval_data_cv(dat = dat, dat1 = dat1, k.cv = k.cv, seeed = seeed, time_inval = c(91.5,1860.5,60.5), ncl = ncl)
    
  }
  if(data_src == "other"){
    
    if(is.null(dat1)) dat1 = data_input
    
    out = eval_data_cv(dat = data_input, dat1 = dat1, k.cv = k.cv, seeed = seeed, time_inval = tint, ncl = ncl)
    
  }

  return(out)
}
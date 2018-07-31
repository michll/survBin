#' Evaluation of the data input using a standard survival method (LASSO) and the classifier version (GBM)
#'
#' @param dat data input for LASSO penalization method
#' @param dat1 data input for classifier approach
#' @param k.cv number of cross-validation folds (default: 10)
#' @param seeed should the seed be fixed (default: 123)
#' @param time_inval time sequence for estimating the time dependent AUC curve, defined as c(start, stop, increments)
#' @param ncl number of cores used for parallelization
#'
#' @return predictins per method
#' @export
#'
#' @examples
eval_data_cv = function(dat, dat1 = NULL, k.cv = 10, seeed = 123, time_inval = c(183,915,30.5), ncl = 1){
  
  #sc = sc2 = sc3 = NULL
  
  data.out = data.out1 = res = list()
  set.seed(seeed)
  fold.outer = crossvalFolds(dat$y[,2],k=k.cv)

  sim_res = parallel::mclapply(1:k.cv, function(xx, ll) {
    
    data.out = data.out1 = ll$dat
    
    data.out$x = ll$dat$x[ll$fold.outer!=xx,]
    data.out$y = ll$dat$y[ll$fold.outer!=xx,]
    
    if(!is.null(ll$dat1)){
      data.out1$x = ll$dat1$x[ll$fold.outer!=xx,]
      data.out1$y = ll$dat1$y[ll$fold.outer!=xx,]
    }
    else{
      data.out1$x = data.out$x
      data.out1$y = data.out$y
    }
    # #--------------
    # # survival glmnet train
    model = glmnet::cv.glmnet(x = as.matrix(data.out$x),y = survival::Surv(data.out$y[,1],data.out$y[,2]),family="cox",alpha=1, nfolds = ll$k.cv)
    # # predict survival
    pred_l = predict(model,newx=as.matrix(ll$dat$x[ll$fold.outer==xx,]),type="link",s="lambda.1se")
    # # time_roc survival
    sc = score_time(time = ll$dat$y[ll$fold.outer==xx,1],death = ll$dat$y[ll$fold.outer==xx,2], pred_l, time_int = ll$time_inval)

        #--------------
    # survival own
    model = surv_fct(x = data.out1$x, y = data.out1$y, tms = seq(as.numeric(ll$time_inval[1]),as.numeric(ll$time_inval[2]),
                                                                 by = as.numeric(ll$time_inval[3])) )
    #data_model = time_rs_fct(tsq = model$tseq, m = model$pred, y = data.out1$y)
    #s1 = survival::coxph(survival::Surv(tstart, tstop, event_new) ~ rs_t , data = data_model)
    # predict survival
    pred_s = predict_surv(model, as.matrix(ll$dat1$x[ll$fold.outer==xx,]))
    #data_test = time_rs_fct(tsq = model$tseq, m = pred_s$pred, y = matrix(NA,ncol = 2, nrow = nrow(dat1$x[fold.outer==i,])))
    #p1 = predict(s1, newdata = data_test[!is.na(data_test$rs_t),], collaps = data_test$id[!is.na(data_test$rs_t)])
    #p1 = p1[order(match(names(p1), paste0("id",1:length(p1))))]
    # time_roc survival
    #sc2 = score_time(time = dat1$y[fold.outer==i,1],death = dat1$y[fold.outer==i,2], p1, time_int = time_inval)
    sc3 = score_time(time = ll$dat1$y[ll$fold.outer==xx,1],death = ll$dat1$y[ll$fold.outer==xx,2], pred_s$pred_int, time_int = ll$time_inval)
    #------------
    # Random survival forest
    datrf = data.frame(time = data.out1$y[,1], event = data.out1$y[,2], data.out1$x)
    model = ranger::ranger( data = datrf, dependent.variable.name = "time", status.variable.name = "event" )
    pred_r = predict(model, ll$dat1$x[ll$fold.outer==xx,])
    pred_r = 1-apply(pred_r$survival,1,median)
    
    sc_rf = score_time(time = ll$dat1$y[ll$fold.outer==xx,1],death = ll$dat1$y[ll$fold.outer==xx,2], pred_r, time_int = ll$time_inval)
    #-------------
    return( c(sc, sc3, sc_rf) )
    #sim_res[[i]] = c(sc, sc3)
    
  }, ll = list(dat = dat, dat1 = dat1, fold.outer = fold.outer, k.cv = k.cv, time_inval = time_inval), mc.cores = ncl)
#  parallel::stopCluster(cl)
  
  out = matrix(unlist(sim_res), nrow = length(sim_res), byrow = T)
  res$pred_lasso = out[,1]
  #res$pred_idsurv = out[,2]
  res$pred_int = out[,2]
  res$pred_rf = out[,3]
  return(res)
  
}


#' Survival modelling as a classification model
#'
#' @param x predictor matrix
#' @param y response as a Surv() object
#' @param classif classification method for caret package (default: gbm)
#' @param tms sequence of time intervals to be evaluated
#'
#' @return sequence of classification models per time point
#' @export
#'
#' @examples
surv_fct = function(x, y, classif = "gbm", tms = NULL){
  if(is.null(tms)){
    tt = c(round(quantile(y[,1],0.1),1),round(quantile(y[,1],0.95),1))
    tseq = seq(tt[1],tt[2],length.out = 20)
  }
  else{
    tseq = tms
  }
  
  datx = x
  m1 = var_imp = res = cv_loss = list()
  pred = matrix(NA, nrow = nrow(datx), ncol = length(tseq))
  for(i in 1:length(tseq)){
    daty = ifelse(y[,1]>tseq[i], 0, ifelse(y[,2]==1, 1,NA))
    if(any(daty == 1, na.rm = T) & any(daty == 0, na.rm = T) & !all(is.na(daty))){
      ii = !is.na(daty)
      if(classif == "gbm"){
        # parametersGrid = expand.grid(
        #   shrinkage = 0.1,
        #   n.minobsinnode = 10,
        #   n.trees=c(100,500,1000),
        #   interaction.depth=c(5, 10, 15)
        # )
        fitControl = caret::trainControl(method="cv", number = 10,classProbs = T, summaryFunction=twoClassSummary, allowParallel = F)
        m1[[i]] = caret::train(x = datx[ii,], y = paste0("X",daty[ii]), method = "gbm",metric = "ROC", 
                               trControl=fitControl, verbose=F)#, tuneGrid=parametersGrid)
        var_imp[[i]] = caret::varImp(m1[[i]])
        cv_loss[[i]] = mean(m1[[i]]$resample$ROC, na.rm = T)
        pred[,i] = predict(m1[[i]], newdata = as.matrix(datx), type = "prob")[,2]
      }
      if(classif == "xgboost"){
        parametersGrid = expand.grid(
          eta = 0.1,
          colsample_bytree=c(0.55, 0.7, 0.85),
          max_depth=3:5,
          nrounds=200,
          gamma=1,
          min_child_weight=c(1, 5, 10),
          subsample=c(0.7, 0.8, 0.9)
        )
        control <- caret::trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=mnLogLoss, allowParallel=F)
        m1[[i]] <- caret::train(x = as.matrix(datx[ii,]), y = paste0("X",daty[ii]), method="xgbTree", metric="logLoss", 
                                trControl=control, tuneGrid=parametersGrid)
        # var_imp[[i]] = caret::varImp(m1[[i]]) # hier xgboost verwenden
        pred[,i] = predict(m1[[i]], newdata = as.matrix(datx), type = "prob")[,2]
      }
      # datx = cbind(datx, pred[,1:i])
    }
    else{
      m1[[i]] = NA
      var_imp[[i]] = NA
      pred[,i] = rep(0,nrow(datx))
    }
    #datx = cbind(datx, pred[,1:i])
  }
  res$cv_loss = cv_loss
  res$var_imp = var_imp
  res$tseq = tseq
  res$model = m1
  res$pred = pred
  #res$surv = coxph(Surv(y[,1],y[,2]) ~ pred)
  return(res)
}
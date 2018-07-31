feat_det = function(model, ncut = 50){

  res = list()
  tsq = model$tseq
  m1 = model$model
  vi = model$var_imp
  vi_mat = matrix(0, ncol = length(tsq), nrow = nrow(vi[[1]]$importance))
  colnames(vi_mat) = paste0("t_",tsq)
  rownames(vi_mat) = rownames(vi[[1]]$importance)
  for(i in 1:length(tsq)){
    vi_mat[,i] = vi[[i]]$importance$Overall
  }
  ii = which(apply(vi_mat,1,max, na.rm = T)>ncut)
  res$tseq = tsq
  res$imp = sort(apply(vi_mat,1,mean, na.rm=T), decreasing = T)
  res$namesii = rownames(vi_mat)[ii]
  res$vi_mat = vi_mat
  res$ind = ii
  class(res) = "feat_det"
  return(res)
}

# dat = imp1
# data_input = dat[dat$STUDYID == "EFC6546",]
# data_input = dat
# not_to_use <- c("RPT", "STUDYID", "LKADT_P", "DEATH", "DISCONT", "ENDTRS_C", "ENTRT_PC")
# not_to_use_Idx <- which(colnames(data_input) %in% not_to_use)
# 
# x <- model.matrix( ~ 0 + ., (data_input[, -not_to_use_Idx])[, apply(data_input[, -not_to_use_Idx], 2, function(l) !any(is.na(l)))])
# 
# dat = res = list()
# dat$y = survival::Surv(time = data_input$LKADT_P, event = data_input$DEATH)
# dat$x = x
# dat1 = dat
# dat1$x = data_input[, -not_to_use_Idx]
# # remove columns with complete NAs
# dat1$x = dat1$x[,apply(dat1$x, 2, function(x) !anyNA(x))]

# model1 = surv_fct(x = dat1$x, y = dat1$y, tms = seq(183,915,by = 30.5))
# model2 = surv_fct(x = dat1$x, y = dat1$y, tms = seq(183,915,by = 30.5))
# 
# # Random survival forest
# datrf = data.frame(time = dat1$y[,1], event = dat1$y[,2], dat1$x)
# model = ranger::ranger( Surv(time ,event) ~ ., data = datrf )
# pred_r = predict(model, as.matrix(ll$dat$x[ll$fold.outer==xx,]))
# pred_r = apply(pred_r$survival,1,median)
# sc2 = score_time(time = ll$dat1$y[ll$fold.outer==xx,1],death = ll$dat1$y[ll$fold.outer==xx,2], pred_r, time_int = ll$time_inval)
# 
# featm1 = feat_det(model1, ncut = 75)
# plot(featm1)
# featm2 = feat_det(model2, ncut = 50)
# debugonce(plot.feat_det)
# plot(featm2)
# 
# col1 = c(wesanderson::wes_palette("Darjeeling"), wesanderson::wes_palette("Moonrise2"))
# featm1$vi_mat[featm1$ind,]
# library(reshape2)

# plot(featm1)



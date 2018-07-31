##################
# Simulate survival data - NEW
##################

library("survival")
library("timeROC")
library("Bolstad2")
library('glmnet')
library('sampling')
library('ipred')
library('doSNOW')
library('foreach')
library('parallel') 
#library('doMC')
library('caret')
library('doParallel')
library('R.utils')
library('xgboost')
R.utils::sourceDirectory("src")





n = c(200, 500)
p = c(50, 1000)
corel = c(0, 0.35, 0.7)
k.cv = 10
rept = 5

ll = expand.grid(n, p, corel)
ll_rep = NULL
for(i in 1:nrow(ll)){
  ll_rep = rbind(ll_rep, matrix(rep(ll[i,],rept), ncol = 3, byrow = T))
}

res = list()
sc = sc1 = sc2 = sc3 = matrix(NA, nrow = k.cv, ncol = nrow(ll_rep))
# set.seed(123)
ncl = 3
# Start cluster as nodes on the kernels
# cl=parallel::makeCluster(ncl,type="SOCK")
# doSNOW::registerDoSNOW(cl)
# Linux
doParallel::registerDoParallel(ncl)
sim_res = foreach::foreach(j = 1:nrow(ll_rep),.errorhandling = "pass") %dopar% {
# sim_res = foreach::foreach(j = 1:3,.errorhandling = "pass") %dopar% {
# for(j in 1:1){
  pj = as.numeric(ll_rep[j,2]); nj = as.numeric(ll_rep[j,1]); corj = as.numeric(ll_rep[j,3])
  nn = pj+nj+corj
  set.seed(nn+j)
  dat = sim_surv(nj, c(1,0.5,-1,-0.5), pj, correl = corj, c(0.5,0.55), ia = 1)
  
  data.out = list()
  fold.outer = crossvalFolds(dat$y[,2],k=k.cv) 
  for(i in 1:k.cv){
    data.out$x = dat$x[fold.outer!=i,]
    data.out$y = dat$y[fold.outer!=i,]
    
    # survival glmnet train
    model = cv.glmnet(x = as.matrix(data.out$x),y = Surv(data.out$y[,1],data.out$y[,2]),family="cox",alpha=1, nfolds = k.cv)
    # predict survival
    pred_l = predict(model,newx=as.matrix(dat$x[fold.outer==i,]),type="link",s="lambda.1se")
    # time_roc survival
    sc[i,j] = score_time(time = dat$y[fold.outer==i,1],death = dat$y[fold.outer==i,2], pred_l, time_int = c(round(quantile(dat$y[,1],0.1),1),round(quantile(dat$y[,1],0.9),1),0.1))
    
    # survival own
    model = surv_fct(x = data.out$x, y = data.out$y)
    data_model = time_rs_fct(tsq = model$tseq, m = model$pred, y = data.out$y)
    s1 = coxph(Surv(tstart, tstop, event_new) ~ rs_t , data = data_model)
    # predict survival
    pred_s = predict_surv(model, as.matrix(dat$x[fold.outer==i,]))
    data_test = time_rs_fct(tsq = model$tseq, m = pred_s$pred, y = matrix(NA,ncol = 2, nrow = nrow(dat$x[fold.outer==i,])))
    p1 = predict(s1, newdata = data_test, collaps = data_test$id)
    p1 = p1[order(match(names(p1), paste0("id",1:length(p1))))]
    # time_roc survival
    sc2[i,j] = score_time(time = dat$y[fold.outer==i,1],death = dat$y[fold.outer==i,2], p1, time_int = c(round(quantile(dat$y[,1],0.1),1),round(quantile(dat$y[,1],0.9),1),0.1))
    sc3[i,j] = score_time(time = dat$y[fold.outer==i,1],death = dat$y[fold.outer==i,2], pred_s$pred_int, time_int = c(round(quantile(dat$y[,1],0.1),1),round(quantile(dat$y[,1],0.9),1),0.1))
    #sc3[i,j] = score_time(time = dat$y[fold.outer==i,1],death = dat$y[fold.outer==i,2], apply(apply(pred_s$pred, 1, weighting), 2, mean), time_int = c(round(quantile(dat$y[,1],0.1),1),round(quantile(dat$y[,1],0.9),1),0.1))
    
  }
  res[[1]] = sc; res[[2]] = sc2; res[[3]] = sc3
}
# parallel::stopCluster(cl)

save(sim_res, file = "sim_res.rda")

vi=NULL
for(i in 1:length(model$var_imp)){
  vi = cbind(vi,model$var_imp[[i]]$importance$Overall)
}
apply(vi,1,mean)


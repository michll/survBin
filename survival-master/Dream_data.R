#########################
# Dream data
#########################

# load data, imp1-imp5 and impLBFI
load(file  ="C:/Users/laimighm/Downloads/Prostate/final/data_new/goldstandard_1a.Rdata")

dat = imp1

levels(dat[, "AGEGRP2"]) <- paste0("cat", 1:3)
levels(dat[, "RACE_C"]) <- paste0("cat", 1:4)    
levels(dat[, "REGION_C"]) <- paste0("cat", 1:5)

# study = dat[dat$STUDYID == "EFC6546",]
# study = dat[dat$STUDYID == "CELGENE",]
# study = dat[dat$STUDYID == "ASCENT2",]

debugonce(dream_data_eval)
celgene_pred = dream_data_eval(data_input = dat[dat$STUDYID == "CELGENE",])
ascent2_pred = dream_data_eval(data_input = dat[dat$STUDYID == "ASCENT2",])
efc6546_pred = dream_data_eval(data_input = dat[dat$STUDYID == "EFC6546",])

cel_asc_pred = dream_data_eval(data_input = dat[dat$STUDYID %in% c("CELGENE","ASCENT2"),])
cel_efc_pred = dream_data_eval(data_input = dat[dat$STUDYID %in% c("CELGENE","EFC6546"),])
asc_efc_pred = dream_data_eval(data_input = dat[dat$STUDYID %in% c("ASCENT2","EFC6546"),])



#################################

# target variables and IDs
not_to_use <- c("RPT", "STUDYID", "LKADT_P", "DEATH", "DISCONT", "ENDTRS_C", "ENTRT_PC")
not_to_use_Idx <- which(colnames(study) %in% not_to_use)

x <- model.matrix( ~ 0 + ., (study[, -not_to_use_Idx])[, apply(study[, -not_to_use_Idx], 2, function(l) !any(is.na(l)))])
dim(x)

dat = list()
dat$y = Surv(time = study$LKADT_P, event = study$DEATH)
dat$x = x
dat1 = dat
dat1$x = study[, -not_to_use_Idx]
# remove columns with complete NAs
dat1$x = dat1$x[,apply(dat1$x, 2, function(x) !anyNA(x))]
k.cv = 10
sc = sc2 = sc3 = NULL

data.out = data.out1 = list()
set.seed(123)
fold.outer = crossvalFolds(dat$y[,2],k=k.cv) 
for(i in 1:k.cv){
  data.out$x = dat$x[fold.outer!=i,]
  data.out$y = dat$y[fold.outer!=i,]
  
  data.out1$x = dat1$x[fold.outer!=i,]
  data.out1$y = dat1$y[fold.outer!=i,]
  
  # survival glmnet train
  model = cv.glmnet(x = as.matrix(data.out$x),y = Surv(data.out$y[,1],data.out$y[,2]),family="cox",alpha=1, nfolds = k.cv)
  # predict survival
  pred_l = predict(model,newx=as.matrix(dat$x[fold.outer==i,]),type="link",s="lambda.1se")
  # time_roc survival
  sc[i] = score_time(time = dat$y[fold.outer==i,1],death = dat$y[fold.outer==i,2], pred_l, time_int = c(183,915,30.5))
  
  # survival own
  #model = surv_fct(x = data.out$x, y = data.out$y, tms = seq(122,915,30.5))
  model = surv_fct(x = data.out1$x, y = data.out1$y, tms = seq(122,915,30.5))
  data_model = time_rs_fct(tsq = model$tseq, m = model$pred, y = data.out1$y)
  s1 = coxph(Surv(tstart, tstop, event_new) ~ rs_t , data = data_model)
  # predict survival
  pred_s = predict_surv(model, as.matrix(dat1$x[fold.outer==i,]))
  data_test = time_rs_fct(tsq = model$tseq, m = pred_s$pred, y = matrix(NA,ncol = 2, nrow = nrow(dat1$x[fold.outer==i,])))
  p1 = predict(s1, newdata = data_test, collaps = data_test$id)
  p1 = p1[order(match(names(p1), paste0("id",1:length(p1))))]
  # time_roc survival
  sc2[i] = score_time(time = dat1$y[fold.outer==i,1],death = dat1$y[fold.outer==i,2], p1, time_int = c(183,915,30.5))
  sc3[i] = score_time(time = dat1$y[fold.outer==i,1],death = dat1$y[fold.outer==i,2], pred_s$pred_int, time_int = c(183,915,30.5))
  #sc3[i,j] = score_time(time = dat$y[fold.outer==i,1],death = dat$y[fold.outer==i,2], apply(apply(pred_s$pred, 1, weighting), 2, mean), time_int = c(round(quantile(dat$y[,1],0.1),1),round(quantile(dat$y[,1],0.9),1),0.1))
  
}


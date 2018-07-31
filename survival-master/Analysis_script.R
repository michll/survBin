###################
# Analysis script
###################

library("survival")
library("timeROC")
library("Bolstad2")
library('glmnet')
library('sampling')
library('ipred')
library('ggplot2')
library('foreach')
library('parallel') 
#library('doMC')
library('caret')
library('doParallel')
#library('R.utils')
library('xgboost')
library('gbm')
library('ranger')
R.utils::sourceDirectory("src")

###################
# Sim data
###################

sim_output3 = data_sel_eval(data_input, k.cv = 10, seeed = 123, data_src = "sim", data_spec = list(n = c(200, 500),
                                                                                     p = c(50, 1000),
                                                                                     corel = c(0, 0.35, 0.7),
                                                                                     beta_vec = c(1,0.5,-1,-0.5),
                                                                                     beta_ia = 1,
                                                                                     surv_cens = c(0.5,0.55),
                                                                                     rept = 5), 
              ncl = 10, tint = NULL)
save(sim_output3, file = "sim_output3.rda")

sim_output = data_sel_eval(data_input, k.cv = 3, seeed = 123, data_src = "sim", data_spec = list(n = c(200),
                                                                                                  p = c(50),
                                                                                                  corel = c(0),
                                                                                                  beta_vec = c(1,0.5,-1,-0.5),
                                                                                                  beta_ia = 1,
                                                                                                  surv_cens = c(0.5,0.55),
                                                                                                  rept = 2), 
                           ncl = 3, tint = NULL)

date()
###################
# Dream data
###################
load(file  ="C:/Users/laimighm/Downloads/Prostate/final/data_new/goldstandard_1a.Rdata")

dat = imp1

levels(dat[, "AGEGRP2"]) <- paste0("cat", 1:3)
levels(dat[, "RACE_C"]) <- paste0("cat", 1:4)    
levels(dat[, "REGION_C"]) <- paste0("cat", 1:5)

celgene_pred = data_sel_eval(data_input = dat[dat$STUDYID == "CELGENE",], data_src = "dream", ncl = 10)
save(celgene_pred, file = "celgene_pred.rda")
ascent2_pred = data_sel_eval(data_input = dat[dat$STUDYID == "ASCENT2",], data_src = "dream", ncl = 10)
save(ascent2_pred, file = "ascent2_pred.rda")
efc6546_pred = data_sel_eval(data_input = dat[dat$STUDYID == "EFC6546",], data_src = "dream", ncl = 10)
save(efc6546_pred, file = "efc6546_pred.rda")

cel_asc_pred = data_sel_eval(data_input = dat[dat$STUDYID %in% c("CELGENE","ASCENT2"),], data_src = "dream", ncl = 10)
save(cel_asc_pred, file = "cel_asc_pred.rda")
cel_efc_pred = data_sel_eval(data_input = dat[dat$STUDYID %in% c("CELGENE","EFC6546"),], data_src = "dream", ncl = 10)
save(cel_efc_pred, file = "cel_efc_pred.rda")
asc_efc_pred = data_sel_eval(data_input = dat[dat$STUDYID %in% c("ASCENT2","EFC6546"),], data_src = "dream", ncl = 10)
save(asc_efc_pred, file = "asc_efc_pred.rda")

all_pred = data_sel_eval(data_input = dat[dat$STUDYID %in% c("ASCENT2","EFC6546","CELGENE"),], data_src = "dream", ncl = 10)
save(all_pred, file = "all_pred.rda")

###################
# TCGA data
###################
TCGA_mutations <- dnet::dRDataLoader(RData='TCGA_mutations')
debugonce(eval_data_cv)
tcga_pred = data_sel_eval(data_input = TCGA_mutations, data_src = "tcga", ncl = 10)

####################
# Gene expr data
####################
source(file = "../../Gene Expr Data/Gene_expr_source_data.R")

date()
pred1456 = data_sel_eval(te1456, data_src = "other", ncl=10, tint = c(3,8,0.5))
save(pred1456, file = "pred1456er.rda")
pred2034 = data_sel_eval(tr2034, data_src = "other", ncl=10, tint = c(3,9,0.5))
date()
save(pred2034, file = "pred2034er.rda")
date()
pred7390 = data_sel_eval(te7390, data_src = "other", ncl=10, tint = c(3,14,0.5))
save(pred7390, file = "pred7390er.rda")
date()

colnames(x) = paste0("C",1:ncol(x))
dat = list(x = x, y = y)
f_par(dat, ncl=1)
debugonce(eval_data_cv)
data_sel_eval(data_input = dat, data_src = "other", ncl = 10, tint = c(0.5,1,0.25), dat1=dat)

######################
# Figure outline
######################
dat = imp1
data_input = dat[dat$STUDYID == "EFC6546",]
data_input = dat
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

set.seed(123)
model1 = surv_fct(x = dat1$x, y = dat1$y, tms = seq(183,915,by = 30.5))
set.seed(123)
model2 = surv_fct(x = te1456$x, y = te1456$y, tms = seq(3,8,by = 0.5))

featm1 = feat_det(model1, ncut = 75)
plot(featm1, scale_x = 30.5*12)
featm2 = feat_det(model1, ncut = 50)
plot(featm2, scale_x = 30.5*12)

featm3 = feat_det(model2, ncut = 75)
plot(featm3)



col1 = c(wesanderson::wes_palette("Darjeeling"), wesanderson::wes_palette("Moonrise2"))
featm1$vi_mat[featm1$ind,]
library(reshape2)

plot(featm1)



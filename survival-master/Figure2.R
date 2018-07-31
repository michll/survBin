##################
# Figure 2
##################
# performance on simulated data

load(file = "Results/sim_output3.rda")
outo = unlist(lapply(sim_output3, function(x){
  x1 = mean(x$pred_lasso)
  x2 = mean(x$pred_rf)
  x3 = mean(x$pred_int)
  
  c(x1,x2,x3)
}))
ll = expand.grid(c(200,500), c(50,1000), c(0,0.35,0.7))
ll_rep = NULL
for(i in 1:nrow(ll)){
  ll_rep = rbind(ll_rep, matrix(rep(ll[i,],5*3), ncol = 3, byrow = T))
}
colnames(ll_rep) = c("N", "p", "correl")

plot_fig1 = data.frame(tAUC = outo, model = factor(rep(c("L1","RF", "BC"), length(sim_output3))), ll_rep)
plot_fig1$N = factor(paste0("N = ",plot_fig1$N), levels = c("N = 200", "N = 500"))
plot_fig1$p = factor(paste0("p = ",plot_fig1$p), levels = c("p = 50", "p = 1000"))
plot_fig1$Np = as.factor(paste0(plot_fig1$N,plot_fig1$p))
plot_fig1$correl = as.factor(paste0("cor = ",plot_fig1$correl))

bp <- ggplot(plot_fig1, aes(x=model, y=tAUC)) + theme_dark(base_size = 12) + labs(x = "") +
  geom_boxplot()  + geom_jitter(position=position_jitter(0.2), col = alpha("red", alpha = 0.25)) + facet_grid(correl~N*p) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red")
bp

ggsave(bp, filename = "Figure_2.pdf", width = 8, height = 8)
##################
# Figure 3
##################
# performance on real data
load(file = "Results/ascent2_pred.rda");
load(file = "Results/celgene_pred.rda");
load(file = "Results/efc6546_pred.rda");
load(file = "Results/asc_efc_pred.rda")
load(file = "Results/cel_asc_pred.rda")
load(file = "Results/cel_efc_pred.rda")
load(file = "Results/all_pred.rda")
outo = c(unlist(ascent2_pred),unlist(celgene_pred),unlist(efc6546_pred), unlist(asc_efc_pred), unlist(cel_asc_pred), unlist(cel_efc_pred), unlist(all_pred) )

plot_fig3 = data.frame(tAUC = outo, model = factor(rep(c(rep("L1",10),rep("BC",10),rep("RF",10)),7),levels = c("L1", "RF","BC")), 
                       study = factor(c(rep("ascent2",30),rep("celgene",30),rep("efc6546",30), rep("asc_efc",30),rep("cel_asc",30),rep("cel_efc",30), rep("all", 30)), 
                                      levels = c("ascent2","celgene","efc6546", "cel_asc", "asc_efc", "cel_efc", "all")) )
bp <- ggplot(plot_fig3, aes(x=model, y=tAUC)) + theme_dark(base_size = 12) + labs(x = "") +
  geom_boxplot()  + geom_jitter(position=position_jitter(0.2), col = alpha("red", alpha = 0.25)) + facet_grid(.~study) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red")
bp

ggsave(bp, filename = "Figure_3A.pdf", width = 8, height = 8)

#################
#Figure 3B
#################
load(file = "Results/pred1456.rda")
load(file = "Results/pred2034.rda")
load(file = "Results/pred7390.rda")
library(ggplot2)

outo = c(unlist(pred1456),unlist(pred2034),unlist(pred7390) )

plot_fig3 = data.frame(tAUC = outo, model = factor(rep(c(rep("L1",10),rep("BC",10),rep("RF",10)),3),levels = c("L1", "RF","BC")), 
                       study = factor(c(rep("s_1456",30),rep("s_2034",30),rep("s_7390",30) ), 
                                      levels = c("s_1456","s_2034","s_7390")) )
bp <- ggplot(plot_fig3, aes(x=model, y=tAUC)) + theme_dark(base_size = 12) + labs(x = "") +
  geom_boxplot()  + geom_jitter(position=position_jitter(0.2), col = alpha("red", alpha = 0.25)) + facet_grid(.~study) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red") + coord_cartesian(ylim=c(0.4, 0.95))
bp

ggsave(bp, filename = "Figure_3B.pdf", width = 8, height = 8)
##################
# Figure 4
##################
# var importance + time dependency + time dep of risk score (model$pred)





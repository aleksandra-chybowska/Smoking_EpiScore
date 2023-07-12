library(ggplot2)
library(ggpubr)
load('/Volumes/marioni-lab/Ola/Lab/Smoking/Trejo_weights/posteriorSummarySM.RData')
pips = as.data.frame(sm_result$PIP)
pips$cpg = rownames(pips)
pips = as.data.frame(pips[-c(1:2),])
colnames(pips) = c("value", "cpg")
return = data.frame("cpg" = pips$cpg, "value" = pips$value)
write.csv(return, "/Volumes/marioni-lab/Ola/Lab/Smoking/Trejo_weights/posteriorSummarySM.csv", row.names = F)

################

load("/Users/shirin/Downloads/meanModelSM.RData")
betas = as.data.frame(meanBetas)
betas = data.frame(cpg = rownames(betas), value = betas$meanBetas)
dupa = as.data.frame(meanBetas)
write.csv(betas, "/Volumes/marioni-lab/Ola/Lab/Smoking/Trejo_weights/meanModelSM.csv", row.names = F)

View(subset(betas, cpg %in% c("age", "Male", "cg05575921", "cg23079012", "cg11207515", "rs2853383_G", "rs2721010_G")))

###############
# Compare Trejo weights with model weights
###############

betas_trejo = read.csv('/Volumes/marioni-lab/Ola/Lab/Smoking/Trejo_weights/meanModelSM.csv')
#betas_ola = read.csv('/Users/shirin/Desktop/pack_years_17865_complete_meanbeta_pip.tsv', sep='\t') # complete
#betas_ola = read.csv('/Users/shirin/Desktop/pack_years_17833_meanbeta_pip.tsv', sep='\t')
betas_ola = read.csv('/Users/shirin/Desktop/pack_years_17833_meanbeta_pip.tsv', sep='\t')

betas_ola = betas_ola[c("Name", "Mean_Beta")]
colnames(betas_ola) = c("cpg", "value_ola")

means = merge(betas_ola, betas_trejo, by="cpg")
means = subset(means, !(value == 0 & value_ola == 0))
ggplot(means, aes(x=value, y=value_ola)) +
  geom_point() +
  geom_smooth(method='lm', colour = "red") +
  xlab("beta Ola") +
  ylab("beta Trejo") +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 7, label.x = -0.2)

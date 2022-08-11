# Read in sweep data with recombination rate

taig_recom_clr <- read.table('taig.ld_recom_converted.w200k_s200k.sweep_sig_olap_counts.bed')
colnames(taig_recom_clr) <- c('scaff','start','stop','taig_cmmb5','taig_cmmb1','taig_cmmb200','sweep_ovlap')
taig_recom_clr$sweep_pres_taig <- ifelse(taig_recom_clr$sweep_ovlap>0,1,0)

coll_recom_clr <- read.table('coll.ld_recom_converted.w200k_s200k.sweep_sig_olap_counts.bed')
colnames(coll_recom_clr) <- c('scaff','start','stop','coll_cmmb5','coll_cmmb1','coll_cmmb200','coll_sweep_ovlap')
coll_recom_clr$sweep_pres_coll <- ifelse(coll_recom_clr$coll_sweep_ovlap>0,1,0)
# filter collared outliers
coll_recom_clr <- coll_recom_clr[coll_recom_clr$coll_cmmb5<20,]

# Load packages
library(plyr)
library(ggplot2)
library(performance)
# Combine collared and taiga data
coll_taig_recom_clr <- join(taig_recom_clr,coll_recom_clr, by = c('scaff','start'), type = 'inner')
coll_taig_recom_clr <- coll_taig_recom_clr[,-9]
head(coll_taig_recom_clr)

## logistic regression for taiga sweeps
taig_taig <- glm(coll_taig_recom_clr$sweep_pres_taig~coll_taig_recom_clr$taig_cmmb5, family="binomial")
taig_coll <- glm(coll_taig_recom_clr$sweep_pres_taig~coll_taig_recom_clr$coll_cmmb5, family="binomial")

## calculate r2
r2_mcfadden(taig_coll)
r2_mcfadden(taig_taig)

## logistic regression for collared sweeps
coll_coll <- glm(coll_taig_recom_clr$sweep_pres_coll~coll_taig_recom_clr$coll_cmmb5, family="binomial")
coll_taig <- glm(coll_taig_recom_clr$sweep_pres_coll~coll_taig_recom_clr$taig_cmmb5, family="binomial")

## calculate r2
r2_mcfadden(coll_coll)
r2_mcfadden(coll_taig)

## create figures for logistic regressions
ggplot(coll_taig_recom_clr, aes(x=taig_cmmb5, y=sweep_pres_taig)) +
  geom_point(color='grey50', alpha =0.75) +
  stat_smooth(method="glm", se = FALSE, method.args = list(family=binomial), color='darkblue', size = 0.75, lty = 2, fullrange=TRUE) +
  xlab('Taiga recombination rate (cM/Mb)')+
  ylab('Selective sweep presence taiga')+
  scale_y_continuous(breaks = c(0,1), labels = c('No','Yes'), limits=c(0,1))+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggplot(coll_taig_recom_clr, aes(x=taig_cmmb5, y=sweep_pres_coll)) +
  geom_point(color='grey50', alpha =0.75) +
  stat_smooth(method="glm", se = FALSE, method.args = list(family=binomial), color='darkblue', size = 0.75, lty = 2, fullrange=TRUE) +
  xlab('Taiga recombination rate (cM/Mb)')+
  ylab('Selective sweep presence collared')+
  scale_y_continuous(breaks = c(0,1), labels = c('No','Yes'), limits=c(0,1))+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


ggplot(coll_taig_recom_clr, aes(x=coll_cmmb5, y=sweep_pres_coll)) +
  geom_point(color='grey50', alpha =0.75) +
  stat_smooth(method="glm", se = FALSE, method.args = list(family=binomial), color='darkblue', size = 0.75, lty = 2, fullrange=TRUE) +
  xlab('Collared recombination rate (cM/Mb)')+
  ylab('Selective sweep presence collared')+
  scale_y_continuous(breaks = c(0,1), labels = c('No','Yes'), limits=c(0,1))+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggplot(coll_taig_recom_clr, aes(x=coll_cmmb5, y=sweep_pres_taig)) +
  geom_point(color='grey50', alpha =0.75) +
  stat_smooth(method="glm", se = FALSE, method.args = list(family=binomial), color='darkblue', size = 0.75, lty = 2, fullrange=TRUE) +
  xlab('Collared recombination rate (cM/Mb)')+
  ylab('Selective sweep presence taiga')+
  scale_y_continuous(breaks = c(0,1), labels = c('No','Yes'), limits=c(0,1))+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


#### Check fst relationship ####

## load fst peak data
fst_peak_class <- read.table('fst_wins.w200k.s200k.all_non_z_scaff.fst_peak_class.bed')
colnames(fst_peak_class) <- c('scaff','start','stp', 'peak_class')

## combine with recombination data for both species
coll_taig_recom_clr_fst <- join(coll_taig_recom_clr, fst_peak_class, by = c('scaff', 'start'), type = 'inner')
head(coll_taig_recom_clr_fst)
coll_taig_recom_clr_fst$peak_class <- factor(x = coll_taig_recom_clr_fst$peak_class, levels = c('No_peak', 'shared_peak', 'PT_unique','CP_unique'))

## create boxplots
ggplot(coll_taig_recom_clr_fst, aes(x=peak_class, y=taig_cmmb5)) +
  geom_boxplot(fill='grey50')+
  xlab('Fst peak category')+
  ylab('Taiga recombination rate (cM/Mb)')+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.ticks.y = element_line(),
        axis.ticks.x = element_line())

ggplot(coll_taig_recom_clr_fst, aes(x=peak_class, y=coll_cmmb5)) +
  geom_boxplot(fill='grey50')+
  xlab('Fst peak category')+
  ylab('Collared recombination rate (cM/Mb)')+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.ticks.y = element_line(),
        axis.ticks.x = element_line())


#### Permutation test ####

# fst peak categories
cats = c('No_peak', 'shared_peak','PT_unique', 'CP_unique')

# get observed mean differences for taiga
obs_diffs = c()
for(i in cats){
  mean_i = mean(coll_taig_recom_clr_fst$taig_cmmb5[coll_taig_recom_clr_fst$peak_class==i])
  print(mean_i)
  for( j in cats){
    mean_j = mean(coll_taig_recom_clr_fst$taig_cmmb5[coll_taig_recom_clr_fst$peak_class==j])
    obs_diffs = append(obs_diffs, abs(mean_i - mean_j))
  }
  
}
obs_diffs

# run permutations
n_perms = 1000

# dataframe to save permuted differences
perm_diffs = data.frame(matrix(vector(), 0, 16,
                       dimnames=list(c(), c("11", "12", "13","14","21","22","23","24","31","32","33","34","41","42","43","44"))),
                stringsAsFactors=F)

# get differences between peak classes for shuffled categories
for(x in 1:n_perms){
  coll_taig_recom_clr_fst$shuff_class <- sample(coll_taig_recom_clr_fst$peak_class, size = length(coll_taig_recom_clr_fst$peak_class))
  perm_diff_set = c()
  for(i in cats){
    mean_i = mean(coll_taig_recom_clr_fst$taig_cmmb5[coll_taig_recom_clr_fst$shuff_class==i])
    for(j in cats){
      
      mean_j =  mean(coll_taig_recom_clr_fst$taig_cmmb5[coll_taig_recom_clr_fst$shuff_class==j])
      perm_diff_set = append(perm_diff_set, abs(mean_i-mean_j))
    }
    
  }
  perm_diffs = rbind(perm_diffs, perm_diff_set)
}

obs_diffs

# get pvalues for each comparison
for(it in 1:16){
  print(length(perm_diffs[,it][perm_diffs[,it] >= obs_diffs[it]])/n_perms)

  
}


# get observed differences in mean recombination rate for collared flycatcher
obs_diffs = c()
for(i in cats){
  mean_i = mean(coll_taig_recom_clr_fst$coll_cmmb5[coll_taig_recom_clr_fst$peak_class==i])
  print(mean_i)
  for( j in cats){
    mean_j = mean(coll_taig_recom_clr_fst$coll_cmmb5[coll_taig_recom_clr_fst$peak_class==j])
    obs_diffs = append(obs_diffs, abs(mean_i - mean_j))
  }
  
}

obs_diffs

# run permutations for collared flycatcher
n_perms = 1000

perm_diffs = data.frame(matrix(vector(), 0, 16,
                               dimnames=list(c(), c("11", "12", "13","14","21","22","23","24","31","32","33","34","41","42","43","44"))),
                        stringsAsFactors=F)

for(x in 1:n_perms){
  coll_taig_recom_clr_fst$shuff_class <- sample(coll_taig_recom_clr_fst$peak_class, size = length(coll_taig_recom_clr_fst$peak_class))
  perm_diff_set = c()
  for(i in cats){
    mean_i = mean(coll_taig_recom_clr_fst$coll_cmmb5[coll_taig_recom_clr_fst$shuff_class==i])
    for(j in cats){
      
      mean_j =  mean(coll_taig_recom_clr_fst$coll_cmmb5[coll_taig_recom_clr_fst$shuff_class==j])
      perm_diff_set = append(perm_diff_set, abs(mean_i-mean_j))
    }
    
  }
  perm_diffs = rbind(perm_diffs, perm_diff_set)
}

obs_diffs

# print pvalues for collared flycatcher group differences
for(it in 1:16){
  print(length(perm_diffs[,it][perm_diffs[,it] >= obs_diffs[it]])/n_perms)
  

}















#### Load data ####

taig_w200_phase1 <- read.table('taig.scaff_combined.phase1.mean_rho.w200k_s200k.txt')
taig_w200_phase2 <- read.table('taig.scaff_combined.phase2.mean_rho.w200k_s200k.txt')

taig_w1Mb_phase1 <- read.table('taig.scaff_combined.phase1.mean_rho.w1Mb_s1MB.txt')
taig_w1Mb_phase2 <- read.table('taig.scaff_combined.phase2.mean_rho.w1Mb_s1MB.txt')

taig_w5Mb_phase1 <- read.table('taig.scaff_combined.phase1.mean_rho.w5Mb.s5Mb.txt')
taig_w5Mb_phase2 <- read.table('taig.scaff_combined.phase2.mean_rho.w5Mb.s5Mb.txt')

linkmap_w200 <- read.table('Chr.Rec.200kb.5kGap.chromPos.txt')
linkmap_w1Mb <- read.table('Chr.Rec.1000kb.5kGap.txt', header = T)
linkmap_w5Mb <- read.table('Chr.Rec.5000kb.5kGap.bed')


#### Load packages ####  
library(plyr)
library(MASS)

#### Filtering 200kb windows taiga ####

## 5030 windows in beginning
head(taig_w200_phase1)
colnames(taig_w200_phase1) <- c('Scaff', 'Start', 'Stop', 'Rho1', 'Wt1', 'SNPs1')
colnames(taig_w200_phase2) <- c('Scaff', 'Start', 'Stop', 'Rho2', 'Wt2', 'SNPs2')

taig_w200_phase1$Size <- (taig_w200_phase1$Stop - taig_w200_phase1$Start)
taig_w200_phase2$Size <- (taig_w200_phase2$Stop - taig_w200_phase2$Start)

length(taig_w200_phase1$Scaff) ## 5030 windows 
sum(taig_w200_phase1$Size) ## 977268559

taig_w200_phase1 <- taig_w200_phase1[taig_w200_phase1$Size==200000,]
taig_w200_phase2 <- taig_w200_phase2[taig_w200_phase2$Size==200000,]

length(taig_w200_phase2$Scaff) ## 4734 windows
sum(taig_w200_phase1$Size) ## 946800000


#plot(taig_w200_phase1$Rho1, taig_w200_phase2$Rho2, las = 1, xlab = 'Taiga rho phase 1', ylab = 'Taiga rho phase 2', col = 'grey50', cex =0.75)

taig_w200_join <- join(taig_w200_phase1, taig_w200_phase2, by = c('Scaff', 'Start'), type ='inner')
taig_w200_join <- taig_w200_join[taig_w200_join$SNPs1>=100,]
length(taig_w200_join$Scaff) ## 4604 windows
sum(taig_w200_join$Size)  ## 920800000

#plot(taig_w200_join$Rho1, taig_w200_join$Rho2, las = 1, xlab = 'Taiga rho phase 1', ylab = 'Taiga rho phase 2', col = 'grey50')

taig_w200_phase_lm <- lm(taig_w200_join$Rho1 ~ taig_w200_join$Rho2)
standard_res <- rstandard(taig_w200_phase_lm)
taig_w200_join$std_res <- standard_res

plot(taig_w200_join$Rho1, taig_w200_join$Rho2, col = ifelse(abs(taig_w200_join$std_res)>9, 'red', 'black'))

cor.test(taig_w200_join$Rho1, taig_w200_join$Rho2) ## 0.97
cor.test(taig_w200_join$Rho1[taig_w200_join$std_res<9], taig_w200_join$Rho2[taig_w200_join$std_res<9]) ## 0.98

taig_w200_join <- taig_w200_join[abs(taig_w200_join$std_res)<=9,]
length(taig_w200_join$Scaff) ##  4602 windows
sum(taig_w200_join$Size) ## 920400000

taig_w200_join$Mean_rho <- rowMeans(cbind(taig_w200_join$Rho1, taig_w200_join$Rho2))

#### Filtering 1Mb windows taiga ####

colnames(taig_w1Mb_phase1) <- c('Scaff', 'Start', 'Stop', 'Rho1', 'Wt1', 'SNPs1')
colnames(taig_w1Mb_phase2) <- c('Scaff', 'Start', 'Stop', 'Rho2', 'Wt2', 'SNPs2')

taig_w1Mb_phase1$Size <- taig_w1Mb_phase1$Stop - taig_w1Mb_phase1$Start
taig_w1Mb_phase2$Size <- taig_w1Mb_phase2$Stop - taig_w1Mb_phase2$Start

length(taig_w1Mb_phase1$Scaff) ## 1133
sum(taig_w1Mb_phase1$Size) ## 977268559


taig_w1Mb_phase1 <- taig_w1Mb_phase1[taig_w1Mb_phase1$Size==1000000,]
taig_w1Mb_phase2 <- taig_w1Mb_phase2[taig_w1Mb_phase2$Size==1000000,]

length(taig_w1Mb_phase1$Scaff) ## 837
sum(taig_w1Mb_phase1$Size) ## 837000000

#plot(taig_w1Mb_phase1$Rho1, taig_w1Mb_phase2$Rho2)

taig_w1Mb_join <- join(taig_w1Mb_phase1, taig_w1Mb_phase2, by = c('Scaff', 'Start'), type = 'inner')
taig_w1Mb_join <- taig_w1Mb_join[taig_w1Mb_join$SNPs1>=500,]
length(taig_w1Mb_join$Scaff) ## 831
sum(taig_w1Mb_join$Size) ## 831000000

#plot(taig_w1Mb_join$Rho1, taig_w1Mb_join$Rho2)


taig_w1Mb_phase_lm <- lm(taig_w1Mb_join$Rho1 ~ taig_w1Mb_join$Rho2)
standard_res <- rstandard(taig_w1Mb_phase_lm)
taig_w1Mb_join$std_res <- standard_res

plot(taig_w1Mb_join$Rho1, taig_w1Mb_join$Rho2, col = ifelse(abs(taig_w1Mb_join$std_res)>9, 'red', 'black'))

cor.test(taig_w1Mb_join$Rho1, taig_w1Mb_join$Rho2) ## 0.97
cor.test(taig_w1Mb_join$Rho1[abs(taig_w1Mb_join$std_res)<9], taig_w1Mb_join$Rho2[abs(taig_w1Mb_join$std_res)<9]) ## 0.99

taig_w1Mb_join <- taig_w1Mb_join[abs(taig_w1Mb_join$std_res)<=9,]
length(taig_w1Mb_join$Scaff) ## 829
sum(taig_w1Mb_join$Size) ## 829000000

taig_w1Mb_join$Mean_rho <- rowMeans(cbind(taig_w1Mb_join$Rho1, taig_w1Mb_join$Rho2))

#### Filtering 5Mb windows taiga ####


colnames(taig_w5Mb_phase1) <- c('Scaff', 'Start', 'Stop', 'Rho1', 'Wt1', 'SNPs1')
colnames(taig_w5Mb_phase2) <- c('Scaff', 'Start', 'Stop', 'Rho2', 'Wt2', 'SNPs2')

taig_w5Mb_phase1$Size <- taig_w5Mb_phase1$Stop - taig_w5Mb_phase1$Start
taig_w5Mb_phase2$Size <- taig_w5Mb_phase2$Stop - taig_w5Mb_phase2$Start
length(taig_w5Mb_phase1$Scaff) ## 394
sum(taig_w5Mb_phase1$Size) ## 977268559

taig_w5Mb_phase1 <- taig_w5Mb_phase1[taig_w5Mb_phase1$Size==5000000,]
taig_w5Mb_phase2 <- taig_w5Mb_phase2[taig_w5Mb_phase2$Size==5000000,]
length(taig_w5Mb_phase1$Scaff) ## 98
sum(taig_w5Mb_phase1$Size) ## 490000000

taig_w5Mb_join <- join(taig_w5Mb_phase1, taig_w5Mb_phase2, by = c('Scaff', 'Start'), type = 'inner')
plot(taig_w5Mb_join$Rho1, taig_w5Mb_join$Rho2)

length(taig_w5Mb_join$Scaff[taig_w5Mb_join$SNPs1<2500]) ## No windows with less than 0.0005 SNP density

taig_w5Mb_phase_lm <- lm(taig_w5Mb_join$Rho1 ~ taig_w5Mb_join$Rho2)
standard_res <- rstandard(taig_w5Mb_phase_lm)
taig_w5Mb_join$std_res <- standard_res

plot(taig_w5Mb_join$Rho1, taig_w5Mb_join$Rho2,  col = ifelse(abs(taig_w5Mb_join$std_res)>9, 'red', 'black'))

taig_w5Mb_join <- taig_w5Mb_join[abs(taig_w5Mb_join$std_res)<=9,]
length(taig_w5Mb_join$Scaff) ## 97
sum(taig_w5Mb_join$Size) ## 485000000

taig_w5Mb_join$Mean_rho <- rowMeans(cbind(taig_w5Mb_join$Rho1, taig_w5Mb_join$Rho2))

#### Estimating Ne with linkage map ####

head(linkmap_w200)
colnames(linkmap_w200) <- c('Chr', 'ChrStart', 'ChrStop', 'Chr2', 'Scaff', 'Start', 'Stop', 'cmMb')

head(linkmap_w1Mb)
colnames(linkmap_w1Mb) <- c('Scaff', 'Id','Chr', 'cmMb','ChrStart','ChrStop','Scaff2','Start','Stop','Ndata')

head(linkmap_w5Mb)
colnames(linkmap_w5Mb) <- c('Scaff','Start','Stop','cmMb')

taig_w200_join$Start <- taig_w200_join$Start+1
taig_w200_recom_join <- join(linkmap_w200, taig_w200_join, by = c('Scaff', 'Start'), type = 'inner')
head(taig_w200_recom_join)

taig_w1Mb_join$Start <- taig_w1Mb_join$Start + 1
taig_w1Mb_recom_join <- join(linkmap_w1Mb, taig_w1Mb_join, by = c('Scaff','Start'), type = 'inner')
head(taig_w1Mb_recom_join)

taig_w5Mb_recom_join <- join(linkmap_w5Mb, taig_w5Mb_join, by = c('Scaff', 'Start'), type = 'inner')
taig_w5Mb_recom_join$cmMb <- as.numeric(taig_w5Mb_recom_join$cmMb)
head(taig_w5Mb_recom_join)

## Plotting relationships

## 200kb windows
plot(taig_w200_recom_join$Mean_rho ~ taig_w200_recom_join$cmMb + 0)
abline(rlm(taig_w200_recom_join$Mean_rho ~ taig_w200_recom_join$cmMb + 0), lty = 2, col = 'blue')
taig_w200_ld_cm <- rlm(taig_w200_recom_join$Mean_rho ~ taig_w200_recom_join$cmMb + 0)
summary(taig_w200_ld_cm)

###  Coefficients:
###  taig_w200_recom_join$cmMb 
###  0.01271239 

cor.test(taig_w200_recom_join$Mean_rho,taig_w200_recom_join$cmMb) ## 0.46

# Ne:
0.01271239/(4*10^-8) # 317809.8


## 1 Mb windows 
plot(taig_w1Mb_recom_join$Mean_rho ~ taig_w1Mb_recom_join$cmMb + 0)
abline(rlm(taig_w1Mb_recom_join$Mean_rho ~ taig_w1Mb_recom_join$cmMb + 0), lty = 2, col = 'blue')
summary(rlm(taig_w1Mb_recom_join$Mean_rho ~ taig_w1Mb_recom_join$cmMb + 0))

#   Coefficients:
####                         Value   Std. Error t value
##taig_w1Mb_recom_join$cmMb  0.0151  0.0003    51.4254

cor.test(taig_w1Mb_recom_join$Mean_rho,taig_w1Mb_recom_join$cmMb) ## 0.55

# Ne:
0.0151/(4*10^-8) ## 377500

## 5Mb windows

plot(taig_w5Mb_recom_join$Mean_rho ~ taig_w5Mb_recom_join$cmMb + 0)
abline(rlm(taig_w5Mb_recom_join$Mean_rho ~ taig_w5Mb_recom_join$cmMb + 0), lty = 2, col = 'blue')
summary(rlm(taig_w5Mb_recom_join$Mean_rho ~ taig_w5Mb_recom_join$cmMb + 0))

## Coefficients:
##                           Value   Std. Error t value
## taig_w5Mb_recom_join$cmMb  0.0183  0.0007    27.9521

cor.test(taig_w5Mb_recom_join$Mean_rho,taig_w5Mb_recom_join$cmMb) ## 0.77

# Ne:
0.0183/(4*10^-8) ## 457500


#### Load collared data #### 

coll_w200_phase1 <- read.table('coll.scaff_combined.phase1.mean_rho.txt')
coll_w200_phase2  <- read.table('coll.scaff_combined.phase2.mean_rho.txt')

coll_w1Mb_phase1 <- read.table('coll.scaff_combined.phase1.mean_rho.w1Mb_s1MB.txt')
coll_w1Mb_phase2  <- read.table('coll.scaff_combined.phase2.mean_rho.w1Mb_s1MB.txt')

coll_w5Mb_phase1 <- read.table('coll.scaff_combined.phase1.mean_rho.w5Mb.s5Mb.txt') 
coll_w5Mb_phase2 <- read.table('coll.scaff_combined.phase2.mean_rho.w5Mb.s5Mb.txt') 


#### Filtering 200kb windows coll ####

## 5035 windows in beginning
head(coll_w200_phase1)
colnames(coll_w200_phase1) <- c('Scaff', 'Start', 'Stop', 'Rho1', 'Wt1', 'SNPs1')
colnames(coll_w200_phase2) <- c('Scaff', 'Start', 'Stop', 'Rho2', 'Wt2', 'SNPs2')


coll_w200_phase1$Size <- (coll_w200_phase1$Stop - coll_w200_phase1$Start)
coll_w200_phase2$Size <- (coll_w200_phase2$Stop - coll_w200_phase2$Start)
length(coll_w200_phase1$Scaff) ## 5035 windows 
sum(coll_w200_phase1$Size) ## 977998999

coll_w200_phase1 <- coll_w200_phase1[coll_w200_phase1$Size==200000,]
coll_w200_phase2 <- coll_w200_phase2[coll_w200_phase2$Size==200000,]
length(coll_w200_phase2$Scaff) ## 4737 windows
sum(coll_w200_phase1$Size) ## 947400000


#plot(coll_w200_phase1$Rho1, coll_w200_phase2$Rho2)

coll_w200_join <- join(coll_w200_phase1, coll_w200_phase2, by = c('Scaff', 'Start'), type ='inner')
coll_w200_join <- coll_w200_join[coll_w200_join$SNPs1>=100,]
length(coll_w200_join$Scaff) ## 4611 windows
sum(coll_w200_join$Size)  ## 922200000

coll_w200_phase_lm <- lm(coll_w200_join$Rho1 ~ coll_w200_join$Rho2)
standard_res <- rstandard(coll_w200_phase_lm)
coll_w200_join$std_res <- standard_res

plot(coll_w200_join$Rho1, coll_w200_join$Rho2, col = ifelse(abs(coll_w200_join$std_res)>9, 'red', 'black'))

cor.test(coll_w200_join$Rho1, coll_w200_join$Rho2) ## 0.98
cor.test(coll_w200_join$Rho1[coll_w200_join$std_res<9], coll_w200_join$Rho2[coll_w200_join$std_res<9]) ## 0.95

coll_w200_join <- coll_w200_join[abs(coll_w200_join$std_res)<=9,]
length(coll_w200_join$Scaff) ##  4603 windows
sum(coll_w200_join$Size) ## 920600000

coll_w200_join$Mean_rho <- rowMeans(cbind(coll_w200_join$Rho1, coll_w200_join$Rho2))
plot(coll_w200_join$Rho1, coll_w200_join$Rho2)

#### Filtering 1Mb windows coll #### 

colnames(coll_w1Mb_phase1) <- c('Scaff', 'Start', 'Stop', 'Rho1', 'Wt1', 'SNPs1')
colnames(coll_w1Mb_phase2) <- c('Scaff', 'Start', 'Stop', 'Rho2', 'Wt2', 'SNPs2')

coll_w1Mb_phase1$Size <- coll_w1Mb_phase1$Stop - coll_w1Mb_phase1$Start
coll_w1Mb_phase2$Size <- coll_w1Mb_phase2$Stop - coll_w1Mb_phase2$Start
length(coll_w1Mb_phase1$Scaff) ## 1135
sum(coll_w1Mb_phase1$Size) ## 977998999


coll_w1Mb_phase1 <- coll_w1Mb_phase1[coll_w1Mb_phase1$Size==1000000,]
coll_w1Mb_phase2 <- coll_w1Mb_phase2[coll_w1Mb_phase2$Size==1000000,]
length(coll_w1Mb_phase1$Scaff) ## 837
sum(coll_w1Mb_phase1$Size) ## 837000000

plot(coll_w1Mb_phase1$Rho1, coll_w1Mb_phase2$Rho2)

coll_w1Mb_join <- join(coll_w1Mb_phase1, coll_w1Mb_phase2, by = c('Scaff', 'Start'), type = 'inner')
coll_w1Mb_join <- coll_w1Mb_join[coll_w1Mb_join$SNPs1>=500,]
length(coll_w1Mb_join$Scaff) ## 832
sum(coll_w1Mb_join$Size) ## 832000000

#plot(coll_w1Mb_join$Rho1, coll_w1Mb_join$Rho2)


coll_w1Mb_phase_lm <- lm(coll_w1Mb_join$Rho1 ~ coll_w1Mb_join$Rho2)
standard_res <- rstandard(coll_w1Mb_phase_lm)
coll_w1Mb_join$std_res <- standard_res

plot(coll_w1Mb_join$Rho1, coll_w1Mb_join$Rho2, col = ifelse(abs(coll_w1Mb_join$std_res)>9, 'red', 'black'))

cor.test(coll_w1Mb_join$Rho1, coll_w1Mb_join$Rho2) ## 0.98
cor.test(coll_w1Mb_join$Rho1[abs(coll_w1Mb_join$std_res)<9], coll_w1Mb_join$Rho2[abs(coll_w1Mb_join$std_res)<9]) ## 0.99

coll_w1Mb_join <- coll_w1Mb_join[abs(coll_w1Mb_join$std_res)<=9,]
length(coll_w1Mb_join$Scaff) ## 829
sum(coll_w1Mb_join$Size) ## 829000000

coll_w1Mb_join$Mean_rho <- rowMeans(cbind(coll_w1Mb_join$Rho1, coll_w1Mb_join$Rho2))


#### Filtering 5Mb windows coll #### 

colnames(coll_w5Mb_phase1) <- c('Scaff', 'Start', 'Stop', 'Rho1', 'Wt1', 'SNPs1')
colnames(coll_w5Mb_phase2) <- c('Scaff', 'Start', 'Stop', 'Rho2', 'Wt2', 'SNPs2')

coll_w5Mb_phase1$Size <- coll_w5Mb_phase1$Stop - coll_w5Mb_phase1$Start
coll_w5Mb_phase2$Size <- coll_w5Mb_phase2$Stop - coll_w5Mb_phase2$Start
length(coll_w5Mb_phase1$Scaff) ## 394
sum(coll_w5Mb_phase1$Size) ## 977268559

coll_w5Mb_phase1 <- coll_w5Mb_phase1[coll_w5Mb_phase1$Size==5000000,]
coll_w5Mb_phase2 <- coll_w5Mb_phase2[coll_w5Mb_phase2$Size==5000000,]
length(coll_w5Mb_phase1$Scaff) ## 98
sum(coll_w5Mb_phase1$Size) ## 490000000

coll_w5Mb_join <- join(coll_w5Mb_phase1, coll_w5Mb_phase2, by = c('Scaff', 'Start'), type = 'inner')
plot(coll_w5Mb_join$Rho1, coll_w5Mb_join$Rho2)

length(coll_w5Mb_join$Scaff[coll_w5Mb_join$SNPs1<2500]) ## No windows with less than 0.0005 SNP density

coll_w5Mb_phase_lm <- lm(coll_w5Mb_join$Rho1 ~ coll_w5Mb_join$Rho2)
standard_res <- rstandard(coll_w5Mb_phase_lm)
coll_w5Mb_join$std_res <- standard_res

plot(coll_w5Mb_join$Rho1, coll_w5Mb_join$Rho2,  col = ifelse(abs(coll_w5Mb_join$std_res)>9, 'red', 'black'))

coll_w5Mb_join <- coll_w5Mb_join[abs(coll_w5Mb_join$std_res)<=9,]
length(coll_w5Mb_join$Scaff) ## 98
sum(coll_w5Mb_join$Size) ## 490000000

coll_w5Mb_join$Mean_rho <- rowMeans(cbind(coll_w5Mb_join$Rho1, coll_w5Mb_join$Rho2))


#### Estimating Ne with linkage map collared ####

head(linkmap_w200)
colnames(linkmap_w200) <- c('Chr', 'ChrStart', 'ChrStop', 'Chr2', 'Scaff', 'Start', 'Stop', 'cmMb')

head(linkmap_w1Mb)
colnames(linkmap_w1Mb) <- c('Scaff', 'Id','Chr', 'cmMb','ChrStart','ChrStop','Scaff2','Start','Stop','Ndata')

head(linkmap_w5Mb)
colnames(linkmap_w5Mb) <- c('Scaff','Start','Stop','cmMb')

coll_w200_join$Start <- coll_w200_join$Start+1
coll_w200_recom_join <- join(linkmap_w200, coll_w200_join, by = c('Scaff', 'Start'), type = 'inner')
head(coll_w200_recom_join)

coll_w1Mb_join$Start <- coll_w1Mb_join$Start + 1
coll_w1Mb_recom_join <- join(linkmap_w1Mb, coll_w1Mb_join, by = c('Scaff','Start'), type = 'inner')
head(coll_w1Mb_recom_join)

coll_w5Mb_recom_join <- join(linkmap_w5Mb, coll_w5Mb_join, by = c('Scaff', 'Start'), type = 'inner')
coll_w5Mb_recom_join$cmMb <- as.numeric(coll_w5Mb_recom_join$cmMb)
head(coll_w5Mb_recom_join)

## Plotting relationships

## 200kb windows
plot(coll_w200_recom_join$Mean_rho ~ coll_w200_recom_join$cmMb + 0)
abline(rlm(coll_w200_recom_join$Mean_rho ~ coll_w200_recom_join$cmMb + 0), lty = 2, col = 'blue')
coll_w200_ld_cm <- rlm(coll_w200_recom_join$Mean_rho ~ coll_w200_recom_join$cmMb + 0)
summary(coll_w200_ld_cm)

### Coefficients:
###                           Value    Std. Error t value 
##coll_w200_recom_join$cmMb   0.0012   0.0000   103.5784

cor.test(coll_w200_recom_join$Mean_rho,coll_w200_recom_join$cmMb) ## 0.48

# Ne:
0.0012/(4*10^-8) # 30000


## 1 Mb windows 
plot(coll_w1Mb_recom_join$Mean_rho ~ coll_w1Mb_recom_join$cmMb + 0)
abline(rlm(coll_w1Mb_recom_join$Mean_rho ~ coll_w1Mb_recom_join$cmMb + 0), lty = 2, col = 'blue')
summary(rlm(coll_w1Mb_recom_join$Mean_rho ~ coll_w1Mb_recom_join$cmMb + 0))

##Coefficients:
###                            Value   Std. Error t value
### coll_w1Mb_recom_join$cmMb  0.0014  0.0000    61.0327

cor.test(coll_w1Mb_recom_join$Mean_rho,coll_w1Mb_recom_join$cmMb) ## 0.53

# Ne:
0.0014/(4*10^-8) ## 35000

## 5Mb windows

plot(coll_w5Mb_recom_join$Mean_rho ~ coll_w5Mb_recom_join$cmMb + 0)
abline(rlm(coll_w5Mb_recom_join$Mean_rho ~ coll_w5Mb_recom_join$cmMb + 0), lty = 2, col = 'blue')
summary(rlm(coll_w5Mb_recom_join$Mean_rho ~ coll_w5Mb_recom_join$cmMb + 0))

##Coefficients:
###                            Value   Std. Error t value
### coll_w5Mb_recom_join$cmMb  0.0017  0.0001    27.8256

cor.test(coll_w5Mb_recom_join$Mean_rho,coll_w5Mb_recom_join$cmMb) ## 0.76

# Ne:
0.0017/(4*10^-8) ## 42500


#### Convert rho taig ####

taig_w200_join$cmMb_5 <- taig_w200_join$Mean_rho/0.0183
taig_w200_join$cmMb_1 <- taig_w200_join$Mean_rho/0.0151
taig_w200_join$cmMb_200 <- taig_w200_join$Mean_rho/0.01271239



#### Convert rho coll ####

coll_w200_join$cmMb_5 <- coll_w200_join$Mean_rho/0.0017
coll_w200_join$cmMb_1 <- coll_w200_join$Mean_rho/0.0014
coll_w200_join$cmMb_200 <- coll_w200_join$Mean_rho/0.0012

#### Subset data ####

taig_w200_recom <- data.frame(as.factor(taig_w200_join$Scaff), taig_w200_join$Start,taig_w200_join$Stop,taig_w200_join$SNPs1, taig_w200_join$Size,taig_w200_join$SNPs2,taig_w200_join$std_res, taig_w200_join$Mean_rho, taig_w200_join$cmMb_5, taig_w200_join$cmMb_1, taig_w200_join$cmMb_200)
colnames(taig_w200_recom) <- c('Scaff','Start','Stop','SNPs1_taig','Size_taig','SNPs2_taig', 'std_res_taig','Mean_rho_taig','cmmb_5_taig','cmmb_1_taig','cmmb_200_taig')

coll_w200_recom <- data.frame(as.factor(coll_w200_join$Scaff),coll_w200_join$Start, coll_w200_join$Stop, coll_w200_join$SNPs1,coll_w200_join$Size,coll_w200_join$SNPs2,coll_w200_join$std_res, coll_w200_join$Mean_rho, coll_w200_join$cmMb_5, coll_w200_join$cmMb_1,coll_w200_join$cmMb_200)
colnames(coll_w200_recom) <- c('Scaff','Start','Stop','SNPs1_coll','Size_coll','SNPs2_coll', 'std_res_coll','Mean_rho_coll','cmmb_5_coll','cmmb_1_coll','cmmb_200_coll')

write.table(taig_w200_recom, 'taig.ld_recom_converted.w200k_s200k.txt', col.names = T, row.names = F, quote=F, sep = '\t')
write.table(coll_w200_recom, 'coll.ld_recom_converted.w200k_s200k.txt', col.names = T, row.names = F, quote=F, sep = '\t')

#### Correlation between species ####

coll_taig_w200_join <- join(coll_w200_recom, taig_w200_recom, by = c('Scaff', 'Start'), type = 'inner')
cor.test(coll_taig_w200_join$cmmb_5_coll[coll_taig_w200_join$cmmb_5_coll<20], coll_taig_w200_join$cmmb_5_taig[coll_taig_w200_join$cmmb_5_coll<20]) ## 0.45
plot(coll_taig_w200_join$cmmb_5_coll, coll_taig_w200_join$cmmb_5_taig)

plot(coll_taig_w200_join$cmmb_5_coll[coll_taig_w200_join$cmmb_5_coll<20], coll_taig_w200_join$cmmb_5_taig[coll_taig_w200_join$cmmb_5_coll<20])






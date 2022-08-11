## Load some packages

library(signal)
library(plyr)
library(dplyr)


## Read the data and perform some filtering
coll_pied_fst_z <- read.table('Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/coll.pied.weir_fst.w200k.s200k.all_z_scaff.windowed.weir.fst.chrompos')
head(coll_pied_fst_z)
colnames(coll_pied_fst_z) <- c('chr','start','stop','sites','wt_fst','fst')
coll_pied_fst_z$mid <- (coll_pied_fst_z$start+coll_pied_fst_z$stop)/2
coll_pied_fst_z <- coll_pied_fst_z[order(coll_pied_fst_z$mid),]
coll_pied_fst_z <- coll_pied_fst_z[coll_pied_fst_z$sites>=200,]
length(coll_pied_fst_z$chr) ##283

par_taig_fst_z <- read.table('Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/par.taig.weir_fst.w200k.s200k.all_z_scaff.windowed.weir.fst.chrompos')
head(par_taig_fst_z)
colnames(par_taig_fst_z)  <- c('chr','start','stop','sites','wt_fst','fst')
par_taig_fst_z$mid <- (par_taig_fst_z$start+par_taig_fst_z$stop)/2
par_taig_fst_z <- par_taig_fst_z[order(par_taig_fst_z$mid),]
par_taig_fst_z <- par_taig_fst_z[par_taig_fst_z$sites>=200,]
length(par_taig_fst_z$chr) ## 284

##### Writing functions necessary to locate peaks

## Calc z_fst function
z_fst <- function(x){
  (x-mean(x))/sd(x)
}

## Function to take zfst and calculate across each chromosome
z_fst_chr <- function(y){
  chr_list <- unique(y[,1])
  chr <- y[,1]
  z_fst_vec <- c()
  for(i in 1:length(chr_list)){
    fst <- y[,5]
    current_chr <- chr_list[i]
    z_fst_chr_vals <- z_fst(fst[chr==current_chr])
    z_fst_vec <- append(z_fst_vec, z_fst_chr_vals)
  }
  return(z_fst_vec)
}

## Write a function to smooth windows by chromosome
## Takes two variables x = values of fst to smooth
## y = dataframe with chromosome values as first column
smooth_windows <- function(x,y){
  chr <- y[,1]
  zfst <- x
  chr_list <- unique(chr)
  
  zfst_smoothed <- c()
  for(m in 1:length(chr_list)){
    current_chr <- chr_list[m]
    zfst_smooth_chr <- sgolayfilt(zfst[chr==current_chr], p = smooth_p, n = smooth_n)
    zfst_smoothed <- append(zfst_smoothed, zfst_smooth_chr)
  }
  return(zfst_smoothed)
}





## Variables that need to be set to run the diff functions:

## Define smoothing function
# p = power
# n = number wins to smooth over
smooth_p <- 3
smooth_n <- 7

## Define number of permutations to perform
n_permutations <- 10000

## Define significance threshold for plotting
p_threshold <- 0.01

## To get the significance of Fst I have a few steps
## 1. Smooth Z-fst across chroms to get observed values
## 2. Run permutations with Z-fst values randomly sorted across genome, recalculate smoothed values
## 3. Compare observed smoothed values to permuted smoothed values to get pvalues
## 4. Run plotting function to check across chroms
## 5. Write output to file


## Get z-fst values for both comps

par_taig_zfst_zchr <- z_fst_chr(par_taig_fst_z)
coll_pied_zfst_zch <- z_fst_chr(coll_pied_fst_z)


## Smooth observed values

par_taig_zfst_smooth <- smooth_windows(par_taig_zfst_zchr, par_taig_fst_z)
coll_pied_zfst_smooth <- smooth_windows(coll_pied_zfst_zch, coll_pied_fst_z)

## Permute windows 

par_taig_zfst_permute <- permute_zfst(par_taig_zfst_zchr, par_taig_fst_z)
coll_pied_zfst_permute <- permute_zfst(coll_pied_zfst_zch, coll_pied_fst_z)

## Get pvalue

par_taig_zfst_pval <- pval_zfst(par_taig_zfst_permute, par_taig_zfst_smooth)
coll_pied_zfst_pval <- pval_zfst(coll_pied_zfst_permute, coll_pied_zfst_smooth)

## Write out significant windows 
#par_taig_sig_wins <- par_taig_fst_z[par_taig_zfst_pval<=p_threshold,]
#coll_pied_sig_wins <- coll_pied_fst_z[coll_pied_zfst_pval<=p_threshold,]

par_taig_sig_wins <- par_taig_fst_z[par_taig_zfst_smooth>2,]
coll_pied_sig_wins <- coll_pied_fst_z[coll_pied_zfst_smooth>2,]

head(par_taig_sig_wins)
head(coll_pied_sig_wins)

length(par_taig_sig_wins$chr) ## 10
length(coll_pied_sig_wins$chr) ## 15

#write.table(par_taig_sig_wins, 'Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/par.taig.weir_fst.w200k.s200k.all_z_scaff.sig_fst_wins.pval_01.txt', col.names = F, row.names = F, quote = F, sep = '\t')
#write.table(coll_pied_sig_wins, 'Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/coll.pied.weir_fst.w200k.s200k.all_z_scaff.sig_fst_wins.pval_01.txt', col.names = F, row.names = F, quote = F, sep = '\t')

write.table(par_taig_sig_wins, 'Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/par.taig.weir_fst.w200k.s200k.all_z_scaff.sig_fst_wins.smooth_zscore_2.txt', col.names = F, row.names = F, quote = F, sep = '\t')
write.table(coll_pied_sig_wins, 'Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/coll.pied.weir_fst.w200k.s200k.all_z_scaff.sig_fst_wins.smooth_zscore_2.txt', col.names = F, row.names = F, quote = F, sep = '\t')

plot(coll_pied_fst_z$mid, sgolayfilt(coll_pied_fst_z$wt_fst, n=smooth_n, p = smooth_p), col = ifelse(coll_pied_zfst_smooth>2, 'red','black'), ylim = c(0,1))
plot(coll_pied_fst_z$mid, sgolayfilt(coll_pied_fst_z$wt_fst, n=smooth_n, p = smooth_p), col = ifelse(coll_pied_zfst_zch>2, 'red','black'), ylim = c(0,1))
plot(coll_pied_fst_z$mid, coll_pied_fst_z$wt_fst, col = ifelse(coll_pied_zfst_smooth>2, 'red','black'), ylim=c(0,1))


plot(par_taig_fst_z$mid, sgolayfilt(par_taig_fst_z$wt_fst, n=smooth_n, p = smooth_p), col = ifelse(par_taig_zfst_pval<=p_threshold, 'red',ifelse(par_taig_zfst_pval<=0.05,'blue','black')), ylim = c(0,1))

#### Autosome data ####

## Read the data and perform some filtering
coll_pied_fst_auto <- read.table('Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/coll.pied.w200k.s200k.all_scaff.non_z_scaff.filtered_sites.windowed.weir.fst.chrompos')
head(coll_pied_fst_auto)
colnames(coll_pied_fst_auto) <- c('chr','start','stop','sites','wt_fst','fst')
coll_pied_fst_auto$mid <- (coll_pied_fst_auto$start+coll_pied_fst_auto$stop)/2
coll_pied_fst_auto <- coll_pied_fst_auto[order(coll_pied_fst_auto$chr,coll_pied_fst_auto$mid),]
coll_pied_fst_auto <- coll_pied_fst_auto[coll_pied_fst_auto$sites>=200,]
head(coll_pied_fst_auto)
length(coll_pied_fst_auto$chr) ## 4589

par_taig_fst_auto <- read.table('Library/Mobile Documents/com~apple~CloudDocs/PhD/research/z_study_chapter_3/fst/par.taig.w200k.s200k.all_scaff.non_z_scaff.filtered_sites.windowed.weir.fst.chrompos')
head(par_taig_fst_auto)
colnames(par_taig_fst_auto)  <- c('chr','start','stop','sites','wt_fst','fst')
par_taig_fst_auto$mid <- (par_taig_fst_auto$start+par_taig_fst_auto$stop)/2
par_taig_fst_auto <- par_taig_fst_auto[order(par_taig_fst_auto$chr, par_taig_fst_auto$mid),]
par_taig_fst_auto <- par_taig_fst_auto[par_taig_fst_auto$sites>=200,]
head(par_taig_fst_auto)
length(par_taig_fst_auto$chr) ## 4629
par_taig_fst_auto <- par_taig_fst_auto[par_taig_fst_auto$chr!='ChrLGE22',]
par_taig_fst_auto <- par_taig_fst_auto[par_taig_fst_auto$chr!='Chr25',]

## Get z-fst values for both comps

par_taig_zfst_auto <- z_fst_chr(par_taig_fst_auto)
coll_pied_zfst_auto <- z_fst_chr(coll_pied_fst_auto)


## Smooth observed values

par_taig_zfst_smooth_auto <- smooth_windows(par_taig_zfst_auto, par_taig_fst_auto)
coll_pied_zfst_smooth_auto <- smooth_windows(coll_pied_zfst_auto, coll_pied_fst_auto)


## Write sig wins

par_taig_sig_wins_auto <- par_taig_fst_auto[par_taig_zfst_smooth_auto>2,]
coll_pied_sig_wins_auto <- coll_pied_fst_auto[coll_pied_zfst_smooth_auto>2,]


write.table(par_taig_sig_wins_auto, 'par.taig.weir_fst.w200k.s200k.all_non_scaff.sig_fst_wins.smooth_zscore_2.txt',col.names = F, row.names = F, quote = F, sep = '\t')
write.table(coll_pied_sig_wins_auto, 'coll.pied.weir_fst.w200k.s200k.all_non_scaff.sig_fst_wins.smooth_zscore_2.txt', col.names = F, row.names = F, quote = F, sep = '\t')











import pandas as pd
import numpy as np
import argparse
from sys import argv

parser = argparse.ArgumentParser(description = 'Calculate pN/pS on a gene by gene basis')
parser.add_argument('-f4', '--four', dest = 'four', help = 'File containing count data for four-fold degenerate sites')
parser.add_argument('-f0' '--zero', dest = 'zero', help = 'File containing count data for zero-fold degenerate sites')
parser.add_argument('-o' '--out', dest = 'out', help = 'Name for output file')

args = parser.parse_args()

zero_fold_file = args.zero
four_fold_file = args.four
out_file = args.out

## Read in count data for zero-fold and four-fold sites
zero_fold = pd.read_csv(zero_fold_file, sep = '\t', names = ['scaff','start','stop','alleles','chrs','allele1','count1','allele2','count2', 'gene'])
four_fold = pd.read_csv(four_fold_file, sep = '\t', names = ['scaff','start','stop','alleles','chrs','allele1','count1','allele2','count2', 'gene'])

# Get allele freq for zero fold sites
zero_fold['p_zero'] = zero_fold['count1']/zero_fold['chrs']
zero_fold['q_zero'] = zero_fold['count2']/zero_fold['chrs']
# get 2pq for each site
zero_fold['het_zero'] = zero_fold['p_zero']*zero_fold['q_zero']*2

# Get allele freq for four fold sites
four_fold['p_four'] = four_fold['count1']/four_fold['chrs']
four_fold['q_four'] = four_fold['count2']/four_fold['chrs']
# Get 2pq for each site
four_fold['het_four'] = four_fold['p_four']*four_fold['q_four']*2

# sum heterozygosity for each gene for zero fold
zero_het_sums = zero_fold.groupby(['gene'])['het_zero'].sum().reset_index(name='het_zero')
# get total sites for each gene for zero-fold
zero_het_sums['len_zero'] = zero_fold.groupby(['gene'])['het_zero'].size().reset_index(name='size')['size']
# get pn for each gene
zero_het_sums['pn'] = zero_het_sums['het_zero']/zero_het_sums['len_zero']

# sum heterozygosity for each gene for four fold
four_het_sums = four_fold.groupby(['gene'])['het_four'].sum().reset_index(name='het_four')
# get total sites for each gene for four-fold
four_het_sums['len_four'] = four_fold.groupby(['gene'])['het_four'].size().reset_index(name='size')['size']
# get ps for each gene
four_het_sums['ps'] = four_het_sums['het_four']/four_het_sums['len_four']

# save pn and ps values to a new df
pn_ps_df = pd.merge(four_het_sums, zero_het_sums, how = 'outer', on = 'gene')

# add pn/ps column to df 
pn_ps_df['pn_ps'] = pn_ps_df['pn']/pn_ps_df['ps']

# replace inf values with nan
pn_ps_df.replace([np.inf, -np.inf], np.nan, inplace=True)

# write values to df
pn_ps_df.to_csv(out_file, sep = '\t',  index = False,na_rep='NULL')
import pandas as pd
import numpy as np
import argparse
from sys import argv

parser = argparse.ArgumentParser(description = 'Normalize GC-conservative dN and dS estimates by proportion of strong and weak sites')
parser.add_argument('-d', '--dn_ds', dest = 'dn_ds', help = 'File with raw dN and dS estimates from bpp for GC-conservative')
parser.add_argument('-o' '--out', dest = 'out', help = 'Output file name')
args = parser.parse_args()


dn_ds_file = args.dn_ds
out_file = args.out

# Read in dN and dS values still separated by strong-strong and weak-weak
dn_ds = pd.read_csv(dn_ds_file, sep = '\t', names = ['gene', 'dn_s_count', 'dn_s_norm','dn_w_count','dn_w_norm', 'ds_s_count','ds_s_norm','ds_w_count','ds_w_norm'])

# estimate weights for dS for S-S and W-W
dn_ds['ds_GC_wt'] = dn_ds['ds_s_norm']/(dn_ds['ds_s_norm']+dn_ds['ds_w_norm'])
dn_ds['ds_AT_wt'] = dn_ds['ds_w_norm']/(dn_ds['ds_s_norm']+dn_ds['ds_w_norm'])

# estimate weights for dN for S-S and W-W
dn_ds['dn_GC_wt'] = dn_ds['dn_s_norm']/(dn_ds['dn_s_norm']+dn_ds['dn_w_norm'])
dn_ds['dn_AT_wt'] = dn_ds['dn_w_norm']/(dn_ds['dn_s_norm']+dn_ds['dn_w_norm'])

# estimate weighted values of dS for S-S and W-W
dn_ds['ds_GC'] = dn_ds['ds_s_count'] * dn_ds['ds_GC_wt']
dn_ds['ds_AT'] = dn_ds['ds_w_count'] * dn_ds['ds_AT_wt']

# estimate weighted values of dN for S-S and W-W
dn_ds['dn_GC'] = dn_ds['dn_s_count'] * dn_ds['dn_GC_wt']
dn_ds['dn_AT'] = dn_ds['dn_w_count'] * dn_ds['dn_AT_wt']

# estimate total dS value for each gene
dn_ds['ds_tot'] = dn_ds['ds_GC'] + dn_ds['ds_AT']

# estimate total dN value for eac gene
dn_ds['dn_tot'] = dn_ds['dn_GC'] + dn_ds['dn_AT']

# estimate total norms for dS
dn_ds['ds_norm_tot'] = dn_ds['ds_s_norm']+dn_ds['ds_w_norm']
dn_ds['dn_norm_tot'] = dn_ds['dn_s_norm']+dn_ds['dn_w_norm']

# estimate dN and dS for each gene, dividing by norms
dn_ds['dn_gene'] = dn_ds['dn_tot']/dn_ds['dn_norm_tot']
dn_ds['ds_gene'] = dn_ds['ds_tot']/dn_ds['ds_norm_tot']

# estimate a gene based dN/dS value
dn_ds['dn_ds_gene'] = dn_ds['dn_gene']/dn_ds['ds_gene']

# in case of any genes having a negative norm estimate from bio++ remove these outliers
dn_ds = dn_ds[dn_ds['dn_s_norm']>=0]
dn_ds = dn_ds[dn_ds['dn_w_norm']>=0]
dn_ds = dn_ds[dn_ds['ds_s_norm']>=0]
dn_ds = dn_ds[dn_ds['ds_w_norm']>=0]

# write output to file
dn_ds.to_csv(out_file, sep = '\t', index = False, na_rep='NULL', columns=['gene','ds_tot','dn_tot','ds_norm_tot','dn_norm_tot','dn_gene','ds_gene','dn_ds_gene'])
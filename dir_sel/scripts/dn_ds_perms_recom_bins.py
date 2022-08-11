import pandas as pd
import numpy as np
import argparse
from sys import argv
import os

parser = argparse.ArgumentParser(description = 'Estimate overall pN/pS, dN/dS and optionally alpha and omega a, for a given set of genes, as well as jackknife resampled standard error.')
parser.add_argument('-g','--gene_tot',dest = 'gene_tot', help = 'List of all genes')
parser.add_argument('-d', '--dn_ds', dest = 'dn_ds', help = 'File containing dN and dS estimates for genes.')
parser.add_argument('-p', '--pn_ps', dest = 'pn_ps', help = 'File containing pN and pS estimates for genes.')
parser.add_argument('-f4', dest='freq_file_four', help = 'File with derived allele counts for neutral sites')
parser.add_argument('-f0', dest='freq_file_zero', help = 'File with derived allele counts for selected sites')
parser.add_argument('-n', dest='n_alleles', help = 'Number of alleles sampled', type = int)
parser.add_argument('--dfe', dest = 'dfe_path', help = 'Path to dfe-alpha executable')
parser.add_argument('-pn', '--perm', dest = 'perm_num', type = int)
parser.add_argument('-t2', dest='time2', type = float, help = 't2 parameter for dfe-alpha')
parser.add_argument('-n2', dest='n2', type = int, help = 'n2 parameter for dfe-alpha')
parser.add_argument('-o' '--out', dest = 'out', help = 'Prefix to use for outfiles')
args = parser.parse_args()

freq_file_four = args.freq_file_four 
freq_file_zero = args.freq_file_zero 
n_alleles = args.n_alleles 
dfe_path = args.dfe_path
out_pre = args.out
time2 = args.time2
n2 = args.n2



dn_ds_file = args.dn_ds
pn_ps_file = args.pn_ps
gene_tot_file = args.gene_tot
perms = args.perm_num

dn_ds_df = pd.read_csv(dn_ds_file, sep = '\t')
pn_ps_df = pd.read_csv(pn_ps_file, sep = '\t')
gene_tot = pd.read_csv(gene_tot_file, names = ['gene','bin_type'], sep = '\t')

# write function to calculate dn/ds
def calc_dn_ds(dn_ds_df, gene_list):
    dn_count = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(gene_list)])
    dn_norm = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(gene_list)])
    dn = dn_count/dn_norm

    ds_count = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(gene_list)])
    ds_norm = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(gene_list)])
    ds = ds_count/ds_norm

    dn_ds = dn/ds
    return dn_ds

# write function to calculate pn/ps
def calc_pn_ps(pn_ps_df, gene_list):
    pn = np.sum(pn_ps_df['het_zero'][pn_ps_df['gene'].isin(gene_list)])/np.sum(pn_ps_df['len_zero'][pn_ps_df['gene'].isin(gene_list)])
    ps = np.sum(pn_ps_df['het_four'][pn_ps_df['gene'].isin(gene_list)])/np.sum(pn_ps_df['len_four'][pn_ps_df['gene'].isin(gene_list)])
    pn_ps = pn/ps
    return pn_ps

# calcuate observed dn/ds for all three recom groups
dn_ds_low_obs = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='low']['gene']))
dn_ds_med_obs = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='med']['gene']))
dn_ds_high_obs = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='high']['gene']))

# calculate observed diff in dn/ds for all three groups  
LM_diff_obs_dn_ds = abs(dn_ds_low_obs - dn_ds_med_obs)
MH_diff_obs_dn_ds = abs(dn_ds_high_obs - dn_ds_med_obs)
LH_diff_obs_dn_ds = abs(dn_ds_high_obs - dn_ds_low_obs)

# calcuate observed pn/ps for all three recom groups
pn_ps_low_obs = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='low']['gene']))
pn_ps_med_obs = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='med']['gene']))
pn_ps_high_obs = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='high']['gene']))

# calculate observed diff in pn/ps for all three groups   
LM_diff_obs_pn_ps = abs(pn_ps_low_obs - pn_ps_med_obs)
MH_diff_obs_pn_ps = abs(pn_ps_high_obs - pn_ps_med_obs)
LH_diff_obs_pn_ps = abs(pn_ps_high_obs - pn_ps_low_obs)

# parameters for running dfe-alpha
site0_string = """
data_path_1     /proj/sllstore2017033/nobackup/work/madeline/chapter2/mol_evol/dfe_data_files/data/
sfs_input_file  {0}
est_dfe_results_dir site_0_2_epoch
site_class 0
fold 1
epochs 2
t2_variable 0
t2 {1}
search_n2 0
n2 {2}
"""

site1_string = """
data_path_1     /proj/sllstore2017033/nobackup/work/madeline/chapter2/mol_evol/dfe_data_files/data/
sfs_input_file  {0}
est_dfe_results_dir site_1_2_epoch
site_class 1
fold 1
epochs 2
est_dfe_demography_results_file site_0_2_epoch/est_dfe.out
mean_s_variable 1
mean_s -0.1
beta_variable 1
beta 0.5
"""

alpha_string = """
data_path_1 /proj/sllstore2017033/nobackup/work/madeline/chapter2/mol_evol/dfe_data_files/data/
divergence_file {0}
est_alpha_omega_results_file alph_2_epoch/alpha_omega.txt
est_dfe_results_file site_1_2_epoch/est_dfe.out
neut_egf_file  site_0_2_epoch/neut_egf.out
sel_egf_file site_1_2_epoch/sel_egf.out
do_jukes_cantor 0
remove_poly 0
"""

## If count data for SFS given then run dfe-alpha
if freq_file_four is not None:

    # Read count data for four-fold and zero-fold sites
    freq_counts_four = pd.read_csv(freq_file_four, sep = '\t', names = ['scaff','start','stop','n_chr','n_der','gene'])
    freq_counts_zero = pd.read_csv(freq_file_zero, sep = '\t', names = ['scaff','start','stop','n_chr','n_der','gene'])

    count_vec_four_low = []
    count_vec_zero_low = []
    count_vec_four_med = []
    count_vec_zero_med = []
    count_vec_four_high = []
    count_vec_zero_high = []

    # Get counts for all three groups
    four_counts_sub_low = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['gene']))]
    zero_counts_sub_low = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['gene']))]
    four_counts_sub_med = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['gene']))]
    zero_counts_sub_med = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['gene']))]
    four_counts_sub_high = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['gene']))]
    zero_counts_sub_high = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['gene']))]

    # estimate sfs for all three groups
    for i in range(0,n_alleles+1):
        count_vec_four_low.append(len(four_counts_sub_low[four_counts_sub_low['n_der']==i]))
        count_vec_zero_low.append(len(zero_counts_sub_low[zero_counts_sub_low['n_der']==i]))
        count_vec_four_med.append(len(four_counts_sub_med[four_counts_sub_med['n_der']==i]))
        count_vec_zero_med.append(len(zero_counts_sub_med[zero_counts_sub_med['n_der']==i]))
        count_vec_four_high.append(len(four_counts_sub_high[four_counts_sub_high['n_der']==i]))
        count_vec_zero_high.append(len(zero_counts_sub_high[zero_counts_sub_high['n_der']==i]))
    
    # write sfs files
    sfs_file_low = out_pre + ".low.sfs_for_dfe.txt"
    f = open(sfs_file_low, "w")
    f.write('1')
    f.write('\n')
    f.write(str(n_alleles))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_zero_low))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_four_low))
    f.close()

    sfs_file_med = out_pre + ".med.sfs_for_dfe.txt"
    f = open(sfs_file_med, "w")
    f.write('1')
    f.write('\n')
    f.write(str(n_alleles))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_zero_med))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_four_med))
    f.close()

    sfs_file_high = out_pre + ".high.sfs_for_dfe.txt"
    f = open(sfs_file_high, "w")
    f.write('1')
    f.write('\n')
    f.write(str(n_alleles))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_zero_high))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_four_high))
    f.close()

    # get divergence data for all three groups
    dn_count_low = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['gene']))])
    dn_norm_low = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['gene']))])
    ds_count_low = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['gene']))])
    ds_norm_low = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['gene']))])

    dn_count_med = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['gene']))])
    dn_norm_med = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['gene']))])
    ds_count_med = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['gene']))])
    ds_norm_med = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['gene']))])

    dn_count_high = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['gene']))])
    dn_norm_high = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['gene']))])
    ds_count_high = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['gene']))])
    ds_norm_high = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['gene']))])

    # write divergence data
    sel_div_low = ['1', dn_norm_low, dn_count_low]
    neut_div_low = ['0', ds_norm_low, ds_count_low]

    div_file_low = out_pre + ".low.div_for_alpha.txt"
    f = open(div_file_low, "w")
    f.write(' '.join(str(d) for d in sel_div_low))
    f.write('\n')
    f.write(' '.join(str(d) for d in neut_div_low))
    f.write('\n')
    f.close()

    sel_div_med = ['1', dn_norm_med, dn_count_med]
    neut_div_med = ['0', ds_norm_med, ds_count_med]

    div_file_med = out_pre + ".med.div_for_alpha.txt"
    f = open(div_file_med, "w")
    f.write(' '.join(str(d) for d in sel_div_med))
    f.write('\n')
    f.write(' '.join(str(d) for d in neut_div_med))
    f.write('\n')
    f.close()

    sel_div_high = ['1', dn_norm_high, dn_count_high]
    neut_div_high = ['0', ds_norm_high, ds_count_high]

    div_file_high = out_pre + ".high.div_for_alpha.txt"
    f = open(div_file_high, "w")
    f.write(' '.join(str(d) for d in sel_div_high))
    f.write('\n')
    f.write(' '.join(str(d) for d in neut_div_high))
    f.write('\n')
    f.close()

    # write dfe-alpha control files
    f=open('alpha_dfe_site_0.ctl', "w")
    f.write(site0_string.format(sfs_file_low, time2, n2))
    f.close()

    f=open('alpha_dfe_site_1.ctl', "w")
    f.write(site1_string.format(sfs_file_low))
    f.close()

    f=open('alpha_dfe_alpha.ctl', "w")
    f.write(alpha_string.format(div_file_low))
    f.close()
   
    # run dfe-alpha for all groups
    site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
    site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
    alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

    os.system(site0_cmd)
    os.system(site1_cmd)
    os.system(alpha_cmd)

    rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{}.low.obs_dfe_site_0.txt'.format(out_pre)
    rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{}.low.obs_dfe_site_1.txt'.format(out_pre)
    rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {}.low.obs_alph_omega.txt'.format(out_pre)
    move_div = 'mv {}.low.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
    move_sfs = 'mv {}.low.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

    os.system(rename_site0)
    os.system(rename_site1)
    os.system(rename_alpha)
    os.system(move_div)
    os.system(move_sfs)

    f=open('alpha_dfe_site_0.ctl', "w")
    f.write(site0_string.format(sfs_file_med, time2, n2))
    f.close()

    f=open('alpha_dfe_site_1.ctl', "w")
    f.write(site1_string.format(sfs_file_med))
    f.close()

    f=open('alpha_dfe_alpha.ctl', "w")
    f.write(alpha_string.format(div_file_med))
    f.close()

    site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
    site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
    alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

    os.system(site0_cmd)
    os.system(site1_cmd)
    os.system(alpha_cmd)

    rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{}.med.obs_dfe_site_0.txt'.format(out_pre)
    rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{}.med.obs_dfe_site_1.txt'.format(out_pre)
    rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {}.med.obs_alph_omega.txt'.format(out_pre)
    move_div = 'mv {}.med.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
    move_sfs = 'mv {}.med.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

    os.system(rename_site0)
    os.system(rename_site1)
    os.system(rename_alpha)
    os.system(move_div)
    os.system(move_sfs)

    f=open('alpha_dfe_site_0.ctl', "w")
    f.write(site0_string.format(sfs_file_high, time2, n2))
    f.close()

    f=open('alpha_dfe_site_1.ctl', "w")
    f.write(site1_string.format(sfs_file_high))
    f.close()

    f=open('alpha_dfe_alpha.ctl', "w")
    f.write(alpha_string.format(div_file_high))
    f.close()

    site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
    site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
    alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

    os.system(site0_cmd)
    os.system(site1_cmd)
    os.system(alpha_cmd)

    rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{}.high.obs_dfe_site_0.txt'.format(out_pre)
    rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{}.high.obs_dfe_site_1.txt'.format(out_pre)
    rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {}.high.obs_alph_omega.txt'.format(out_pre)
    move_div = 'mv {}.high.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
    move_sfs = 'mv {}.high.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

    os.system(rename_site0)
    os.system(rename_site1)
    os.system(rename_alpha)
    os.system(move_div)
    os.system(move_sfs)


LM_diff_dn_ds_perm = []
MH_diff_dn_ds_perm  = []
LH_diff_dn_ds_perm  = []

LM_diff_pn_ps_perm = []
MH_diff_pn_ps_perm  = []
LH_diff_pn_ps_perm  = []

LM_diff_alph_perm = []
MH_diff_alph_perm  = []
LH_diff_alph_perm  = []

LM_diff_es_perm = []
MH_diff_es_perm  = []
LH_diff_es_perm  = []

LM_diff_omega_perm = []
MH_diff_omega_perm  = []
LH_diff_omega_perm  = []

# get observed differences for dfe-alpha stats
if freq_file_four is not None:
    # read in observed dfe-alpha out files    
    obs_alpha_file_low = out_pre + '.low.obs_alph_omega.txt'
    obs_dfe_file_low = 'dfe_alpha_out_files/' + out_pre + '.low.obs_dfe_site_1.txt'

    obs_alpha_file_med = out_pre + '.med.obs_alph_omega.txt'
    obs_dfe_file_med = 'dfe_alpha_out_files/' + out_pre + '.med.obs_dfe_site_1.txt'
    obs_alpha_file_high = out_pre + '.high.obs_alph_omega.txt'
    obs_dfe_file_high = 'dfe_alpha_out_files/' + out_pre + '.high.obs_dfe_site_1.txt'


    obs_dfe_low = pd.read_csv(obs_dfe_file_low, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])
    obs_dfe_med = pd.read_csv(obs_dfe_file_med, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])
    obs_dfe_high = pd.read_csv(obs_dfe_file_high, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])

    obs_alph_low = pd.read_csv(obs_alpha_file_low, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
    obs_alph_med = pd.read_csv(obs_alpha_file_med, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
    obs_alph_high = pd.read_csv(obs_alpha_file_high, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])

    # get observed stats
    obs_es_low = obs_dfe_low['es_val'][0]
    obs_es_med = obs_dfe_med['es_val'][0]
    obs_es_high = obs_dfe_high['es_val'][0]
    
    obs_alpha_val_low = obs_alph_low['a2'][0]
    obs_alpha_val_med = obs_alph_med['a2'][0]
    obs_alpha_val_high = obs_alph_high['a2'][0]

    obs_omega_val_low = obs_alph_low['o2'][0]
    obs_omega_val_med = obs_alph_med['o2'][0]
    obs_omega_val_high = obs_alph_high['o2'][0]

    # get observed differences between groups
    LM_diff_es_obs = abs(obs_es_low - obs_es_med)
    MH_diff_es_obs  = abs(obs_es_high - obs_es_med)
    LH_diff_es_obs  = abs(obs_es_high - obs_es_low)

    LM_diff_alpha_obs = abs(obs_alpha_val_low-obs_alpha_val_med)
    MH_diff_alpha_obs  = abs(obs_alpha_val_high-obs_alpha_val_med)
    LH_diff_alpha_obs = abs(obs_alpha_val_high-obs_alpha_val_low)

    LM_diff_omega_obs = abs(obs_omega_val_low-obs_omega_val_med)
    MH_diff_omega_obs = abs(obs_omega_val_high-obs_omega_val_med) 
    LH_diff_omega_obs = abs(obs_omega_val_high-obs_omega_val_low) 

# run permutation tests
for x in range(0, perms):
    # randomly re-shuffle order of gene ids    
    gene_tot['perm_gene'] = np.random.permutation(gene_tot['gene'].values)

    # get re-shuffled values of dn/ds for each group
    low_dn_ds_perm = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))
    med_dn_ds_perm = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))
    high_dn_ds_perm = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))
    
    # get re-shuffled values of pn/ps for each group
    low_pn_ps_perm = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))
    med_pn_ps_perm = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))
    high_pn_ps_perm = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))

    # calculate re-shuffled differences in dn/ds
    LM_diff_dn_ds_perm.append(abs(low_dn_ds_perm - med_dn_ds_perm))
    MH_diff_dn_ds_perm.append(abs(med_dn_ds_perm - high_dn_ds_perm))
    LH_diff_dn_ds_perm.append(abs(high_dn_ds_perm - low_dn_ds_perm))

    # calculate re-shuffled differences in pn/ps
    LM_diff_pn_ps_perm.append(abs(low_pn_ps_perm - med_pn_ps_perm))
    MH_diff_pn_ps_perm.append(abs(med_pn_ps_perm - high_pn_ps_perm))
    LH_diff_pn_ps_perm.append(abs(high_pn_ps_perm - low_pn_ps_perm))

    # run permutations for dfe-lapha
    if freq_file_four is not None:
        count_vec_four_low = []
        count_vec_zero_low = []
        count_vec_four_med = []
        count_vec_zero_med = []
        count_vec_four_high = []
        count_vec_zero_high = []

        # get count data for re-shuffled groups
        four_counts_sub_low = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))]
        zero_counts_sub_low = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))]
        four_counts_sub_med = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))]
        zero_counts_sub_med = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))]
        four_counts_sub_high = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))]
        zero_counts_sub_high = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))]

        # estimate sfs for re-shuffled groups
        for i in range(0,n_alleles+1):
            count_vec_four_low.append(len(four_counts_sub_low[four_counts_sub_low['n_der']==i]))
            count_vec_zero_low.append(len(zero_counts_sub_low[zero_counts_sub_low['n_der']==i]))
            count_vec_four_med.append(len(four_counts_sub_med[four_counts_sub_med['n_der']==i]))
            count_vec_zero_med.append(len(zero_counts_sub_med[zero_counts_sub_med['n_der']==i]))
            count_vec_four_high.append(len(four_counts_sub_high[four_counts_sub_high['n_der']==i]))
            count_vec_zero_high.append(len(zero_counts_sub_high[zero_counts_sub_high['n_der']==i]))
        
        # write re-shuffled sfs to files
        sfs_file_low = out_pre + ".low.sfs_for_dfe.txt"
        f = open(sfs_file_low, "w")
        f.write('1')
        f.write('\n')
        f.write(str(n_alleles))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_zero_low))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_four_low))
        f.close()

        sfs_file_med = out_pre + ".med.sfs_for_dfe.txt"
        f = open(sfs_file_med, "w")
        f.write('1')
        f.write('\n')
        f.write(str(n_alleles))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_zero_med))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_four_med))
        f.close()

        sfs_file_high = out_pre + ".high.sfs_for_dfe.txt"
        f = open(sfs_file_high, "w")
        f.write('1')
        f.write('\n')
        f.write(str(n_alleles))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_zero_high))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_four_high))
        f.close()

        # get reshuffled divergence data
        dn_count_low = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))])
        dn_norm_low = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))])
        ds_count_low = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))])
        ds_norm_low = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='low']['perm_gene']))])

        dn_count_med = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))])
        dn_norm_med = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))])
        ds_count_med = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))])
        ds_norm_med = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='med']['perm_gene']))])

        dn_count_high = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))])
        dn_norm_high = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))])
        ds_count_high = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))])
        ds_norm_high = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='high']['perm_gene']))])


        # write re-shuffled divergence data to files
        sel_div_low = ['1', dn_norm_low, dn_count_low]
        neut_div_low = ['0', ds_norm_low, ds_count_low]

        div_file_low = out_pre + ".low.div_for_alpha.txt"
        f = open(div_file_low, "w")
        f.write(' '.join(str(d) for d in sel_div_low))
        f.write('\n')
        f.write(' '.join(str(d) for d in neut_div_low))
        f.write('\n')
        f.close()

        sel_div_med = ['1', dn_norm_med, dn_count_med]
        neut_div_med = ['0', ds_norm_med, ds_count_med]

        div_file_med = out_pre + ".med.div_for_alpha.txt"
        f = open(div_file_med, "w")
        f.write(' '.join(str(d) for d in sel_div_med))
        f.write('\n')
        f.write(' '.join(str(d) for d in neut_div_med))
        f.write('\n')
        f.close()

        sel_div_high = ['1', dn_norm_high, dn_count_high]
        neut_div_high = ['0', ds_norm_high, ds_count_high]

        div_file_high = out_pre + ".high.div_for_alpha.txt"
        f = open(div_file_high, "w")
        f.write(' '.join(str(d) for d in sel_div_high))
        f.write('\n')
        f.write(' '.join(str(d) for d in neut_div_high))
        f.write('\n')
        f.close()
        
        # write dfe-alpha control files
        f=open('alpha_dfe_site_0.ctl', "w")
        f.write(site0_string.format(sfs_file_low, time2, n2))
        f.close()

        f=open('alpha_dfe_site_1.ctl', "w")
        f.write(site1_string.format(sfs_file_low))
        f.close()

        f=open('alpha_dfe_alpha.ctl', "w")
        f.write(alpha_string.format(div_file_low))
        f.close()
    

        site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
        site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
        alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

        os.system(site0_cmd)
        os.system(site1_cmd)
        os.system(alpha_cmd)

        rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.low.obs_dfe_site_0.txt'.format(out_pre,x)
        rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.low.obs_dfe_site_1.txt'.format(out_pre,x)
        rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {0}.{1}_bs_rep.low.obs_alph_omega.txt'.format(out_pre,x)
        move_div = 'mv {}.low.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
        move_sfs = 'mv {}.low.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

        os.system(rename_site0)
        os.system(rename_site1)
        os.system(rename_alpha)
        os.system(move_div)
        os.system(move_sfs)

        f=open('alpha_dfe_site_0.ctl', "w")
        f.write(site0_string.format(sfs_file_med, time2, n2))
        f.close()

        f=open('alpha_dfe_site_1.ctl', "w")
        f.write(site1_string.format(sfs_file_med))
        f.close()

        f=open('alpha_dfe_alpha.ctl', "w")
        f.write(alpha_string.format(div_file_med))
        f.close()

        site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
        site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
        alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

        os.system(site0_cmd)
        os.system(site1_cmd)
        os.system(alpha_cmd)

        rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.med.obs_dfe_site_0.txt'.format(out_pre,x)
        rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.med.obs_dfe_site_1.txt'.format(out_pre,x)
        rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {0}.{1}_bs_rep.med.obs_alph_omega.txt'.format(out_pre,x)
        move_div = 'mv {}.med.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
        move_sfs = 'mv {}.med.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

        os.system(rename_site0)
        os.system(rename_site1)
        os.system(rename_alpha)
        os.system(move_div)
        os.system(move_sfs)

        f=open('alpha_dfe_site_0.ctl', "w")
        f.write(site0_string.format(sfs_file_high, time2, n2))
        f.close()

        f=open('alpha_dfe_site_1.ctl', "w")
        f.write(site1_string.format(sfs_file_high))
        f.close()

        f=open('alpha_dfe_alpha.ctl', "w")
        f.write(alpha_string.format(div_file_high))
        f.close()

        site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
        site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
        alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

        os.system(site0_cmd)
        os.system(site1_cmd)
        os.system(alpha_cmd)

        rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.high.obs_dfe_site_0.txt'.format(out_pre,x)
        rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.high.obs_dfe_site_1.txt'.format(out_pre,x)
        rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {0}.{1}_bs_rep.high.obs_alph_omega.txt'.format(out_pre,x)
        move_div = 'mv {}.high.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
        move_sfs = 'mv {}.high.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

        os.system(rename_site0)
        os.system(rename_site1)
        os.system(rename_alpha)
        os.system(move_div)
        os.system(move_sfs)

        # read re-shuffled dfe-alpha out files
        boot_alpha_file_low = '{0}.{1}_bs_rep.low.obs_alph_omega.txt'.format(out_pre,x)
        boot_dfe_file_low = 'dfe_alpha_out_files/{0}.{1}_bs_rep.low.obs_dfe_site_1.txt'.format(out_pre,x)

        boot_alpha_file_med = '{0}.{1}_bs_rep.med.obs_alph_omega.txt'.format(out_pre,x)
        boot_dfe_file_med = 'dfe_alpha_out_files/{0}.{1}_bs_rep.med.obs_dfe_site_1.txt'.format(out_pre,x)

        boot_alpha_file_high = '{0}.{1}_bs_rep.high.obs_alph_omega.txt'.format(out_pre,x)
        boot_dfe_file_high = 'dfe_alpha_out_files/{0}.{1}_bs_rep.high.obs_dfe_site_1.txt'.format(out_pre,x)

        boot_alphas_low = pd.read_csv(boot_alpha_file_low, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
        boot_dfes_low = pd.read_csv(boot_dfe_file_low, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])

        boot_alphas_med = pd.read_csv(boot_alpha_file_med, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
        boot_dfes_med = pd.read_csv(boot_dfe_file_med, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])

        boot_alphas_high = pd.read_csv(boot_alpha_file_high, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
        boot_dfes_high = pd.read_csv(boot_dfe_file_high, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])

        # get re-shuffled values
        boot_alph_val_low = boot_alphas_low['a2'][0]
        boot_alph_val_med = boot_alphas_med['a2'][0]
        boot_alph_val_high = boot_alphas_high['a2'][0]

        boot_es_val_low = boot_dfes_low['es_val'][0]
        boot_es_val_med = boot_dfes_med['es_val'][0]
        boot_es_val_high = boot_dfes_high['es_val'][0]

        boot_omega_val_low = boot_alphas_low['o2'][0]
        boot_omega_val_med = boot_alphas_med['o2'][0]
        boot_omega_val_high = boot_alphas_high['o2'][0]

        # calculated re-shuffled differences 
        LM_diff_alph_perm.append(abs(boot_alph_val_low-boot_alph_val_med))
        MH_diff_alph_perm.append(abs(boot_alph_val_high-boot_alph_val_med))
        LH_diff_alph_perm.append(abs(boot_alph_val_high-boot_alph_val_low))

        LM_diff_es_perm.append(abs(boot_es_val_low-boot_es_val_med))
        MH_diff_es_perm.append(abs(boot_es_val_high-boot_es_val_med))
        LH_diff_es_perm.append(abs(boot_es_val_low-boot_es_val_high))

        LM_diff_omega_perm.append(abs(boot_omega_val_low-boot_omega_val_med))
        MH_diff_omega_perm.append(abs(boot_omega_val_high-boot_omega_val_med))
        LH_diff_omega_perm.append(abs(boot_omega_val_high-boot_omega_val_low))


# get p-values for dn/ds for each group comparison
p_dn_ds_LM = sum(i >= LM_diff_obs_dn_ds for i in LM_diff_dn_ds_perm)/perms
p_dn_ds_MH = sum(i >= MH_diff_obs_dn_ds for i in MH_diff_dn_ds_perm)/perms
p_dn_ds_LH = sum(i >= LH_diff_obs_dn_ds for i in LH_diff_dn_ds_perm)/perms

# get p-values for pn/ps for each group comparison
p_pn_ps_LM = sum(i >= LM_diff_obs_pn_ps for i in LM_diff_pn_ps_perm)/perms
p_pn_ps_MH = sum(i >= MH_diff_obs_pn_ps for i in MH_diff_pn_ps_perm)/perms
p_pn_ps_LH = sum(i >= LH_diff_obs_pn_ps for i in LH_diff_pn_ps_perm)/perms

# get p-values for dfe-alpha stats
if freq_file_four is not None:
    p_alph_LM = sum(i >= LM_diff_alpha_obs for i in LM_diff_alph_perm)/perms
    p_alph_MH = sum(i >= MH_diff_alpha_obs for i in MH_diff_alph_perm)/perms
    p_alph_LH = sum(i >= LH_diff_alpha_obs for i in LH_diff_alph_perm)/perms

    p_omega_LM = sum(i >= LM_diff_omega_obs for i in LM_diff_omega_perm)/perms
    p_omega_MH = sum(i >= MH_diff_omega_obs for i in MH_diff_omega_perm)/perms
    p_omega_LH = sum(i >= LH_diff_omega_obs for i in LH_diff_omega_perm)/perms

    p_es_LM = sum(i >= LM_diff_es_obs for i in LM_diff_es_perm)/perms
    p_es_MH = sum(i >= MH_diff_es_obs for i in MH_diff_es_perm)/perms
    p_es_LH = sum(i >= LH_diff_es_obs for i in LH_diff_es_perm)/perms

# write p-values to file
outfile_name = out_pre + '.pvals_out.txt'
if freq_file_four is not None:
    f = open(outfile_name, "w")
    f.write('\t'.join(['Stat','LM','MH','LH']))
    f.write('\n')
    f.write('\t'.join(['dn_ds',str(p_dn_ds_LM),str(p_dn_ds_MH),str(p_dn_ds_LH)]))
    f.write('\n')
    f.write('\t'.join(['pn_ps',str(p_pn_ps_LM),str(p_pn_ps_MH),str(p_pn_ps_LH)]))
    f.write('\n')
    f.write('\t'.join(['alpha', str(p_alph_LM),str(p_alph_MH),str(p_alph_LH)]))
    f.write('\n')
    f.write('\t'.join(['omega', str(p_omega_LM),str(p_omega_MH),str(p_omega_LH)]))
    f.write('\n')
    f.write('\t'.join(['es', str(p_es_LM),str(p_es_MH),str(p_es_LH)]))
    f.close()
else:
    f = open(outfile_name, "w")
    f.write('\t'.join(['Stat','LM','MH','LH']))
    f.write('\n')
    f.write('\t'.join(['dn_ds',str(p_dn_ds_LM),str(p_dn_ds_MH),str(p_dn_ds_LH)]))
    f.write('\n')
    f.write('\t'.join(['pn_ps',str(p_pn_ps_LM),str(p_pn_ps_MH),str(p_pn_ps_LH)]))
    f.close()
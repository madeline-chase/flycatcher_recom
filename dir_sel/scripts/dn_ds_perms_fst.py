import pandas as pd
import numpy as np
import argparse
from sys import argv
import os

parser = argparse.ArgumentParser(description = 'Estimate overall pN/pS, dN/dS and optionally alpha and omega a, for a given set of genes, as well as jackknife resampled standard error.')
parser.add_argument('-g','--gene_tot',dest = 'gene_tot', help = 'List of all genes. Gene name and group ID')
parser.add_argument('-d', '--dn_ds', dest = 'dn_ds', help = 'File containing dN and dS estimates for genes.')
parser.add_argument('-p', '--pn_ps', dest = 'pn_ps', help = 'File containing pN and pS estimates for genes.')
parser.add_argument('-f4', dest='freq_file_four', help = 'File with derived allele counts for neutral sites')
parser.add_argument('-f0', dest='freq_file_zero', help = 'File with derived allele counts for selected sites')
parser.add_argument('-n', dest='n_alleles', help = 'Number of alleles sampled', type = int)
parser.add_argument('--dfe', dest = 'dfe_path', help = 'Path to dfe-alpha executable')
parser.add_argument('-o' '--out', dest = 'out', help = 'Prefix to use for outfiles')
parser.add_argument('-t2', dest='time2', type = float, help = 't2 parameter for dfe-alpha')
parser.add_argument('-n2', dest='n2', type = int, help = 'n2 parameter for dfe-alpha')
parser.add_argument('-pn', '--perm', dest = 'perm_num', type = int)
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

## write function to calculate dn/ds
def calc_dn_ds(dn_ds_df, gene_list):
    dn_count = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(gene_list)])
    dn_norm = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(gene_list)])
    dn = dn_count/dn_norm

    ds_count = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(gene_list)])
    ds_norm = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(gene_list)])
    ds = ds_count/ds_norm

    dn_ds = dn/ds
    return dn_ds

## write function to calculate pn/ps
def calc_pn_ps(pn_ps_df, gene_list):
    pn = np.sum(pn_ps_df['het_zero'][pn_ps_df['gene'].isin(gene_list)])/np.sum(pn_ps_df['len_zero'][pn_ps_df['gene'].isin(gene_list)])
    ps = np.sum(pn_ps_df['het_four'][pn_ps_df['gene'].isin(gene_list)])/np.sum(pn_ps_df['len_four'][pn_ps_df['gene'].isin(gene_list)])
    pn_ps = pn/ps
    return pn_ps

# estimate observed dn/ds for both groups
dn_ds_group1_obs = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='group1']['gene']))
dn_ds_group2_obs = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='group2']['gene']))

# estimate observed difference in dn/ds 
diff_obs_dn_ds = abs(dn_ds_group1_obs - dn_ds_group2_obs)

# estimate observed pn/ps for both groups
pn_ps_group1_obs = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='group1']['gene']))
pn_ps_group2_obs = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='group2']['gene']))

# estimate observed difference in pn/ps   
diff_obs_pn_ps = abs(pn_ps_group1_obs - pn_ps_group2_obs)

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

    count_vec_four_group1 = []
    count_vec_zero_group1 = []
    count_vec_four_group2 = []
    count_vec_zero_group2= []

    ## Subset sites for genes in list
    four_counts_sub_group1 = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['gene']))]
    zero_counts_sub_group1 = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['gene']))]
    four_counts_sub_group2 = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['gene']))]
    zero_counts_sub_group2 = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['gene']))]

    ## Estimate SFS for both groups
    for i in range(0,n_alleles+1):
        count_vec_four_group1.append(len(four_counts_sub_group1[four_counts_sub_group1['n_der']==i]))
        count_vec_zero_group1.append(len(zero_counts_sub_group1[zero_counts_sub_group1['n_der']==i]))
        count_vec_four_group2.append(len(four_counts_sub_group2[four_counts_sub_group2['n_der']==i]))
        count_vec_zero_group2.append(len(zero_counts_sub_group2[zero_counts_sub_group2['n_der']==i]))
    
    ## write SFS files
    sfs_file_group1 = out_pre + ".group1.sfs_for_dfe.txt"
    f = open(sfs_file_group1, "w")
    f.write('1')
    f.write('\n')
    f.write(str(n_alleles))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_zero_group1))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_four_group1))
    f.close()

    sfs_file_group2 = out_pre + ".group2.sfs_for_dfe.txt"
    f = open(sfs_file_group2, "w")
    f.write('1')
    f.write('\n')
    f.write(str(n_alleles))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_zero_group2))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_four_group2))
    f.close()

    # get divergence data
    dn_count_group1 = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['gene']))])
    dn_norm_group1 = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['gene']))])
    ds_count_group1 = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['gene']))])
    ds_norm_group1 = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['gene']))])

    dn_count_group2 = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['gene']))])
    dn_norm_group2 = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['gene']))])
    ds_count_group2 = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['gene']))])
    ds_norm_group2 = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['gene']))])

    # write divergence data files
    sel_div_group1 = ['1', dn_norm_group1, dn_count_group1]
    neut_div_group1 = ['0', ds_norm_group1, ds_count_group1]

    div_file_group1 = out_pre + ".group1.div_for_alpha.txt"
    f = open(div_file_group1, "w")
    f.write(' '.join(str(d) for d in sel_div_group1))
    f.write('\n')
    f.write(' '.join(str(d) for d in neut_div_group1))
    f.write('\n')
    f.close()

    sel_div_group2 = ['1', dn_norm_group2, dn_count_group2]
    neut_div_group2 = ['0', ds_norm_group2, ds_count_group2]

    div_file_group2 = out_pre + ".group2.div_for_alpha.txt"
    f = open(div_file_group2, "w")
    f.write(' '.join(str(d) for d in sel_div_group2))
    f.write('\n')
    f.write(' '.join(str(d) for d in neut_div_group2))
    f.write('\n')
    f.close()

    # write dfe-alpha control files
    f=open('alpha_dfe_site_0.ctl', "w")
    f.write(site0_string.format(sfs_file_group1, time2, n2))
    f.close()

    f=open('alpha_dfe_site_1.ctl', "w")
    f.write(site1_string.format(sfs_file_group1))
    f.close()

    f=open('alpha_dfe_alpha.ctl', "w")
    f.write(alpha_string.format(div_file_group1))
    f.close()
   
    # run dfe-alpha for first group
    site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
    site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
    alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

    os.system(site0_cmd)
    os.system(site1_cmd)
    os.system(alpha_cmd)

    rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{}.group1.obs_dfe_site_0.txt'.format(out_pre)
    rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{}.group1.obs_dfe_site_1.txt'.format(out_pre)
    rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {}.group1.obs_alph_omega.txt'.format(out_pre)
    move_div = 'mv {}.group1.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
    move_sfs = 'mv {}.group1.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

    os.system(rename_site0)
    os.system(rename_site1)
    os.system(rename_alpha)
    os.system(move_div)
    os.system(move_sfs)

    # write dfe-alpha control files for second group
    f=open('alpha_dfe_site_0.ctl', "w")
    f.write(site0_string.format(sfs_file_group2, time2, n2))
    f.close()

    f=open('alpha_dfe_site_1.ctl', "w")
    f.write(site1_string.format(sfs_file_group2))
    f.close()

    f=open('alpha_dfe_alpha.ctl', "w")
    f.write(alpha_string.format(div_file_group2))
    f.close()

    # run dfe-alpha for second group
    site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
    site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
    alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

    os.system(site0_cmd)
    os.system(site1_cmd)
    os.system(alpha_cmd)

    rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{}.group2.obs_dfe_site_0.txt'.format(out_pre)
    rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{}.group2.obs_dfe_site_1.txt'.format(out_pre)
    rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {}.group2.obs_alph_omega.txt'.format(out_pre)
    move_div = 'mv {}.group2.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
    move_sfs = 'mv {}.group2.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

    os.system(rename_site0)
    os.system(rename_site1)
    os.system(rename_alpha)
    os.system(move_div)
    os.system(move_sfs)


diff_dn_ds_perm = []

diff_pn_ps_perm = []

diff_alph_perm = []

diff_es_perm = []

diff_omega_perm = []

# get observed difference for dfe-alpha outputs
if freq_file_four is not None:

    # read dfe-alpha out for group 1
    obs_alpha_file_group1 = out_pre + '.group1.obs_alph_omega.txt'
    obs_dfe_file_group1 = 'dfe_alpha_out_files/' + out_pre + '.group1.obs_dfe_site_1.txt'

    # read dfe-alpha out for group 2
    obs_alpha_file_group2 = out_pre + '.group2.obs_alph_omega.txt'
    obs_dfe_file_group2 = 'dfe_alpha_out_files/' + out_pre + '.group2.obs_dfe_site_1.txt'

    obs_dfe_group1 = pd.read_csv(obs_dfe_file_group1, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])
    obs_dfe_group2 = pd.read_csv(obs_dfe_file_group2, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])

    obs_alph_group1 = pd.read_csv(obs_alpha_file_group1, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
    obs_alph_group2 = pd.read_csv(obs_alpha_file_group2, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])

    # get observed values for both groups
    obs_es_group1 = obs_dfe_group1['es_val'][0]
    obs_es_group2 = obs_dfe_group2['es_val'][0]
    
    obs_alpha_val_group1 = obs_alph_group1['a2'][0]
    obs_alpha_val_group2 = obs_alph_group2['a2'][0]

    obs_omega_val_group1 = obs_alph_group1['o2'][0]
    obs_omega_val_group2 = obs_alph_group2['o2'][0]

    # calculate observed difference between groups
    diff_es_obs = abs(obs_es_group1 - obs_es_group2)

    diff_alpha_obs = abs(obs_alpha_val_group1-obs_alpha_val_group2)

    diff_omega_obs = abs(obs_omega_val_group1-obs_omega_val_group2)

# run permutation tests
for x in range(0, perms):
    # randomly shuffle order of gene ids
    gene_tot['perm_gene'] = np.random.permutation(gene_tot['gene'].values)
    # estimate re-shuffled dn/ds for both groups
    group1_dn_ds_perm = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))
    group2_dn_ds_perm = calc_dn_ds(dn_ds_df, list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))
    
    # estimate re-shuffled pn/ps for both groups
    group1_pn_ps_perm = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))
    group2_pn_ps_perm = calc_pn_ps(pn_ps_df, list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))

    # calculate re-shuffled difference in dn/ds
    diff_dn_ds_perm.append(abs(group1_dn_ds_perm - group2_dn_ds_perm))

    # calculate re-shuffled difference in pn/ps
    diff_pn_ps_perm.append(abs(group1_pn_ps_perm - group2_pn_ps_perm))

    # estimate re-shuffled dfe
    if freq_file_four is not None:
        count_vec_four_group1 = []
        count_vec_zero_group1 = []
        count_vec_four_group2 = []
        count_vec_zero_group2 = []

        # get count data for re-shuffled genes
        four_counts_sub_group1 = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))]
        zero_counts_sub_group1 = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))]
        four_counts_sub_group2 = freq_counts_four[freq_counts_four['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))]
        zero_counts_sub_group2 = freq_counts_zero[freq_counts_zero['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))]

        ## estimate sfs
        for i in range(0,n_alleles+1):
            count_vec_four_group1.append(len(four_counts_sub_group1[four_counts_sub_group1['n_der']==i]))
            count_vec_zero_group1.append(len(zero_counts_sub_group1[zero_counts_sub_group1['n_der']==i]))
            count_vec_four_group2.append(len(four_counts_sub_group2[four_counts_sub_group2['n_der']==i]))
            count_vec_zero_group2.append(len(zero_counts_sub_group2[zero_counts_sub_group2['n_der']==i]))
        
        # write sfs to files
        sfs_file_group1 = out_pre + ".group1.sfs_for_dfe.txt"
        f = open(sfs_file_group1, "w")
        f.write('1')
        f.write('\n')
        f.write(str(n_alleles))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_zero_group1))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_four_group1))
        f.close()

        sfs_file_group2 = out_pre + ".group2.sfs_for_dfe.txt"
        f = open(sfs_file_group2, "w")
        f.write('1')
        f.write('\n')
        f.write(str(n_alleles))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_zero_group2))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_four_group2))
        f.close()

        # write re-shuffled divergence to files
        dn_count_group1 = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))])
        dn_norm_group1 = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))])
        ds_count_group1 = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))])
        ds_norm_group1 = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group1']['perm_gene']))])

        dn_count_group2 = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))])
        dn_norm_group2 = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))])
        ds_count_group2 = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))])
        ds_norm_group2 = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(list(gene_tot[gene_tot['bin_type']=='group2']['perm_gene']))])

        sel_div_group1 = ['1', dn_norm_group1, dn_count_group1]
        neut_div_group1 = ['0', ds_norm_group1, ds_count_group1]

        div_file_group1 = out_pre + ".group1.div_for_alpha.txt"
        f = open(div_file_group1, "w")
        f.write(' '.join(str(d) for d in sel_div_group1))
        f.write('\n')
        f.write(' '.join(str(d) for d in neut_div_group1))
        f.write('\n')
        f.close()

        sel_div_group2 = ['1', dn_norm_group2, dn_count_group2]
        neut_div_group2 = ['0', ds_norm_group2, ds_count_group2]

        div_file_med = out_pre + ".group2.div_for_alpha.txt"
        f = open(div_file_med, "w")
        f.write(' '.join(str(d) for d in sel_div_group2))
        f.write('\n')
        f.write(' '.join(str(d) for d in neut_div_group2))
        f.write('\n')
        f.close()
        
        # write control files and run dfe-alpha for group 1
        f=open('alpha_dfe_site_0.ctl', "w")
        f.write(site0_string.format(sfs_file_group1, time2, n2))
        f.close()

        f=open('alpha_dfe_site_1.ctl', "w")
        f.write(site1_string.format(sfs_file_group1))
        f.close()

        f=open('alpha_dfe_alpha.ctl', "w")
        f.write(alpha_string.format(div_file_group1))
        f.close()
    

        site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
        site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
        alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

        os.system(site0_cmd)
        os.system(site1_cmd)
        os.system(alpha_cmd)

        rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.group1.obs_dfe_site_0.txt'.format(out_pre,x)
        rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.group1.obs_dfe_site_1.txt'.format(out_pre,x)
        rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {0}.{1}_bs_rep.group1.obs_alph_omega.txt'.format(out_pre,x)
        move_div = 'mv {}.group1.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
        move_sfs = 'mv {}.group1.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

        os.system(rename_site0)
        os.system(rename_site1)
        os.system(rename_alpha)
        os.system(move_div)
        os.system(move_sfs)

        # write control files and run dfe-alpha for group 2
        f=open('alpha_dfe_site_0.ctl', "w")
        f.write(site0_string.format(sfs_file_group2, time2, n2))
        f.close()

        f=open('alpha_dfe_site_1.ctl', "w")
        f.write(site1_string.format(sfs_file_group2))
        f.close()

        f=open('alpha_dfe_alpha.ctl', "w")
        f.write(alpha_string.format(div_file_group2))
        f.close()

        site0_cmd = '{0}/est_dfe -c alpha_dfe_site_0.ctl'.format(dfe_path)
        site1_cmd = '{0}/est_dfe -c alpha_dfe_site_1.ctl'.format(dfe_path)
        alpha_cmd = '{0}/est_alpha_omega -c alpha_dfe_alpha.ctl'.format(dfe_path)

        os.system(site0_cmd)
        os.system(site1_cmd)
        os.system(alpha_cmd)

        rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.group2.obs_dfe_site_0.txt'.format(out_pre,x)
        rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{0}.{1}_bs_rep.group2.obs_dfe_site_1.txt'.format(out_pre,x)
        rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {0}.{1}_bs_rep.group2.obs_alph_omega.txt'.format(out_pre,x)
        move_div = 'mv {}.group2.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
        move_sfs = 'mv {}.group2.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

        os.system(rename_site0)
        os.system(rename_site1)
        os.system(rename_alpha)
        os.system(move_div)
        os.system(move_sfs)

        # read dfe-alpha output for re-shuffled data
        boot_alpha_file_group1 = '{0}.{1}_bs_rep.group1.obs_alph_omega.txt'.format(out_pre,x)
        boot_dfe_file_group1 = 'dfe_alpha_out_files/{0}.{1}_bs_rep.group1.obs_dfe_site_1.txt'.format(out_pre,x)

        boot_alpha_file_group2 = '{0}.{1}_bs_rep.group2.obs_alph_omega.txt'.format(out_pre,x)
        boot_dfe_file_group2 = 'dfe_alpha_out_files/{0}.{1}_bs_rep.group2.obs_dfe_site_1.txt'.format(out_pre,x)

        boot_alphas_group1 = pd.read_csv(boot_alpha_file_group1, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
        boot_dfes_group1 = pd.read_csv(boot_dfe_file_group1, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])

        boot_alphas_group2 = pd.read_csv(boot_alpha_file_group2, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
        boot_dfes_group2 = pd.read_csv(boot_dfe_file_group2, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])

        ## get re-shuffled statistics
        boot_alph_val_group1 = boot_alphas_group1['a2'][0]
        boot_alph_val_group2 = boot_alphas_group2['a2'][0]

        boot_es_val_group1 = boot_dfes_group1['es_val'][0]
        boot_es_val_group2= boot_dfes_group2['es_val'][0]

        boot_omega_val_group1 = boot_alphas_group1['o2'][0]
        boot_omega_val_group2 = boot_alphas_group2['o2'][0]

        ## calcuate re-shuffled group differences
        diff_alph_perm.append(abs(boot_alph_val_group1-boot_alph_val_group2))

        diff_es_perm.append(abs(boot_es_val_group1-boot_es_val_group2))

        diff_omega_perm.append(abs(boot_omega_val_group1-boot_omega_val_group2))


# get p-values for dn/ds
p_dn_ds = (sum(i >= diff_obs_dn_ds for i in diff_dn_ds_perm)+1)/perms

# get p-values for pn/ps
p_pn_ps = (sum(i >= diff_obs_pn_ps for i in diff_pn_ps_perm)+1)/perms

# get p-values for dfe-alpha stats
if freq_file_four is not None:
    p_alph= sum(i >= diff_alpha_obs for i in diff_alph_perm)/perms

    p_omega = sum(i >= diff_omega_obs for i in diff_omega_perm)/perms

    p_es = sum(i >= diff_es_obs for i in diff_es_perm)/perms


outfile_name = out_pre + '.pvals_out.txt'

# write out p-values
if freq_file_four is not None:
    f = open(outfile_name, "w")
    f.write('\t'.join(['Stat','12']))
    f.write('\n')
    f.write('\t'.join(['dn_ds',str(p_dn_ds)]))
    f.write('\n')
    f.write('\t'.join(['pn_ps',str(p_pn_ps)]))
    f.write('\n')
    f.write('\t'.join(['alpha', str(p_alph)]))
    f.write('\n')
    f.write('\t'.join(['omega', str(p_omega)]))
    f.write('\n')
    f.write('\t'.join(['es', str(p_es)]))
    f.close()
else:
    f = open(outfile_name, "w")
    f.write('\t'.join(['Stat','12']))
    f.write('\n')
    f.write('\t'.join(['dn_ds',str(p_dn_ds)]))
    f.write('\n')
    f.write('\t'.join(['pn_ps',str(p_pn_ps)]))
    f.close()


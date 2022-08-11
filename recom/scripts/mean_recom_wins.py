import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Write average recombination rate in user defined windows')
parser.add_argument('-r', dest='recom_file', help = 'File with recombination rates of SNP pairs')
parser.add_argument('-w', dest='scaff_file', help = 'File with window coordinates: Scaff Start Stop. 1-based')
parser.add_argument('-s', dest = 'scaffold', help = 'Name of scaffold for analysis')

args = parser.parse_args()

recom_file = args.recom_file
win_file = args.scaff_file
scaffold = args.scaffold

## Read window coordinate file
win_data = pd.read_csv(win_file, sep = '\t', names = ['Scaff', 'Start', 'Stop'])

## Read SNP-pair recombination rate file
recom_data = pd.read_csv(recom_file, sep = ' ', names = ['Left_snp', 'Right_snp', 'Run1', 'Run2','Run3', 'Run4', 'Run5'])

## Take average recombination rate across 5 iterations of rjMCMC
recom_data['Mean_rho'] = recom_data.iloc[:, 2:6].mean(axis=1)

## Convert start to zero base
win_data['Start'] = win_data['Start']-1

## Subset dataset for one scaffold
scaff_win = win_data.copy()
scaff_win = scaff_win.loc[scaff_win.Scaff==scaffold]

## Estimate average rho for each window
for window in range(0, len(scaff_win)):
    win_start = scaff_win.iloc[window]['Start']
    win_stop = scaff_win.iloc[window]['Stop']
    
    ## Get recombination rate data for current window 
    recom_win = recom_data.copy()
    recom_win = recom_win[(recom_win.Right_snp >= win_start) & (recom_win.Left_snp < win_stop)]
    recom_win = recom_win.reset_index(drop=True)

    ## create distance column
    recom_win['dist'] = np.nan
    
    ## Calculate distance between SNP pairs
    ## Distance should be between current SNP and window start if the left SNP is in another window
    ## Distance should be between current SNP and window stop if the right SNP is in another window
    recom_win.loc[recom_win.Left_snp<win_start, 'dist'] = recom_win['Right_snp'] - win_start
    recom_win.loc[recom_win.Right_snp>win_stop, 'dist'] = win_stop - recom_win['Left_snp'] 
    recom_win.loc[(recom_win.Right_snp<=win_stop) & (recom_win.Left_snp>=win_start), 'dist'] = recom_win['Right_snp'] - recom_win['Left_snp']

    ## Get weight for each SNP
    recom_win['weight'] = recom_win['dist']/sum(recom_win['dist'])
    
    ## Calculate weighted mean recombination for window
    window_mean = sum(recom_win['Mean_rho']*recom_win['weight'])
    print(scaffold, win_start, win_stop, window_mean, sum(recom_win['dist']), len(recom_win['dist']), sep = '\t')

import pandas as pd
import argparse
from sys import argv
import os


parser = argparse.ArgumentParser(description = 'Identify sweep regions from adjacent sig values.')
parser.add_argument('-i', '--in', dest = 'input', help = 'SF2 output file')
parser.add_argument('-s', '--scaff', dest = 'scaff', help = 'Name of scaffold')
parser.add_argument('-o', '--out', dest = 'out', help = 'Prefix for outfile')
args = parser.parse_args()

scaff = args.scaff
out_pre = args.out
in_file = args.input

## Read in SF2 output
sf_in = pd.read_csv(in_file, sep = '\t')

## Add column to DF with scaffold for analysis
sf_in.loc[:, 'scaff'] = scaff

## Subset DF for sites with sig CLR values based on sim threshold
sf_sig_df = sf_in[sf_in['LR']>46.25]

## If any sites on the scaffold are significant perform extra filter steps
if len(sf_sig_df)>0:

    ## Write significant sites out to file, saving index from original DF
    sig_df_out = out_pre + ".sig_pos_w_idx.txt"
    sf_sig_df.to_csv(sig_df_out, sep = '\t', index = True, header = False)

    ## Command to merge adjacent significant sites 
    ## Create 'pseudo-positions' using index from DF to merge two significant positions
    ## if they are adjacent in the SF2 output. Report min and max genome coordinate,
    ## total number of sites and mean CLR value. 
    bed_string = "awk -v OFS='\t' '{{print $5,$1,$1+1,$2,$3,$4}}' {0} | bedtools merge -i stdin -c 4,4,4,5 -o min,max,count,mean | awk -v OFS='\t' '{{print $1,$4-1,$5,$6,$7, $5-($4-1)}}' > {1}.sweep_coords_adj_pos_merged.bed".format(sig_df_out, out_pre)
    os.system(bed_string)

    ## Read in file with merged sweep coordinates
    sweep_coord_file = out_pre + ".sweep_coords_adj_pos_merged.bed"
    sweep_coords = pd.read_csv(sweep_coord_file, names = ['scaff','start','stop','count','mean_clr','size'], sep = '\t')

    ## Filter for sweeps with more than one site included and site density greater than
    ## one site/kb
    sweep_coords_sub = sweep_coords[sweep_coords['size']>1]
    sweep_coords_sub =sweep_coords_sub[sweep_coords_sub['count']/sweep_coords_sub['size']>0.001]

    ## Save output for filtered significant sweep coordinates
    sig_df_filt_out = out_pre + ".sweep_coords_adj_pos_merged.filt.bed"
    sweep_coords_sub.to_csv(sig_df_filt_out, index = False, header = True, sep = '\t')
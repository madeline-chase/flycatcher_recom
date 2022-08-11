import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Convert from shapeit hap output to ldhelmet seq file input.')
parser.add_argument('--hap', dest = 'hap_file', help = 'Shapit hap output file')
parser.add_argument('--samp', dest = 'samp_file', help = 'Sample output from shapeit, containing names of samples in the hap file (in order)')
parser.add_argument('--keep', dest = 'keep_file', help = 'File with subset of samples from hap file to keep in the output')

args = parser.parse_args()

hap_file = args.hap_file
samp_file = args.samp_file
keep_file = args.keep_file

# Generate haplotype names for all samples
samp_list = []
for line in open(samp_file):
    samp_name = line.strip().split()[0]
    samp_list.append(samp_name + "_A")
    samp_list.append(samp_name + "_B")

# Remove first four entries from list that don't give sample names
samp_list = samp_list[4:]

# create header to read in phased haplotypes
hap_header = ['Col1', 'Scaff_pos','Pos', 'Allele1','Allele2'] + samp_list

## Create list of samples that should be output
keep_list = []
for line in open(keep_file):
    samp_name = line.strip().split()[0]
    keep_list.append(samp_name + "_A")
    keep_list.append(samp_name + "_B")

# Read in phased data
hap_in = pd.read_csv(hap_file, names = hap_header, sep = ' ')

## Go through phased data to generate sequences
## Add allele 1 nucleotide to seq if 0
## Add allele 2 nucleotide to seq if 1
for samp in keep_list:
    print("> " + samp)
    for pos in range(0,len(hap_in)):
        if hap_in[samp][pos] == 0:
            print(hap_in['Allele1'][pos], end = '', sep='')
        elif hap_in[samp][pos] == 1:
            print(hap_in['Allele2'][pos], end = '', sep='')
    print('\n', end = '')

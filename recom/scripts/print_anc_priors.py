import pandas as pd
from sys import argv

map_file = argv[1]
polar_file = argv[2]

map_df = pd.read_csv(map_file, sep = '\t', names = ['col1', 'scaff_site', 'col3', 'site'])
polar_df = pd.read_csv(polar_file, sep = '\t', names = ['Scaff', 'Site', 'Anc_allele'])

polar_sites = list(polar_df['Site'])

for pos in range(0,len(map_df)):
    current_site = map_df.iloc[pos]['site']
    if current_site in polar_sites:
        anc_allele = polar_df[polar_df['Site']==current_site]['Anc_allele'].item()
        if anc_allele == 'A':
            print(str(current_site -1), '0.97', '0.01', '0.01', '0.01', sep = '\t')
        elif anc_allele == 'C':
            print(str(current_site -1),  '0.01', '0.97','0.01', '0.01', sep = '\t')
        elif anc_allele == 'G':
            print(str(current_site -1),  '0.01', '0.01', '0.97', '0.01', sep = '\t')
        elif anc_allele == 'T':
            print(str(current_site -1),  '0.01', '0.01',  '0.01', '0.97', sep = '\t')
    elif current_site not in polar_sites:
        print(str(current_site-1), '0.25', '0.25', '0.25', '0.25', sep = '\t')
            
            
            
            
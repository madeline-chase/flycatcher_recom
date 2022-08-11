from Bio import SeqIO
from Bio.Seq import Seq
from sys import argv
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Get zero-fold and four-fold degenerate sites from fasta')
parser.add_argument('-f', '--fasta', dest = 'fasta_in', help = 'Fasta formatted file containing coding sequences. CDS must be correct orientation.')
parser.add_argument('-b' '--bed', dest = 'bed_in', help = 'Bed12 formatted file containing coordinates of coding sequences. 0-based.')
parser.add_argument('-o', '--out', dest = 'out', help = 'Outfile name. 1-based positions are written out.')
parser.add_argument('-s', '--scaffs', dest = 'scaffs', help = 'File with chr names, scaffold names size, and orientations. Format: Chr, Scaff, Size, Orient')
args = parser.parse_args()

fasta_file = args.fasta_in
cds_in = args.bed_in
out_file = args.out
scaff_file = args.scaffs

# Read CDS into seq dictionary
fasta_in = open(fasta_file)
record_dict = SeqIO.to_dict(SeqIO.parse(fasta_in, 'fasta'))

# Read cds coordinate file
cds_coords = pd.read_csv(cds_in, sep = '\t', names = ['Scaff','Start','End', 'Name', 'Score', 'Strand', 'ThickStart', 'ThickEnd', 'itemRgb', 'BlockCount', 'BlockSizes', 'BlockStarts'])
## Blocks refer to exons in the CDS
## Block starts is a list of the start position in the gene of each exon, starting at 0 for the start pos of the first exon
## Block sizes gives the length of each exon


# Read scaffold orientation file
scaff_orient = pd.read_csv(scaff_file, sep = '\t', names = ['Chr', 'Scaff', 'Size', 'Orient'])

# Codon table
codon_table = {'TTT': 'Phe', 'TTC': 'Phe','TTA':'Leu','TTG':'Leu','CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu','ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met','GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val', 'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser','CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro','ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr','GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala','TAT':'Tyr','TAC':'Tyr','TAA':'Ter','TAG':'Ter','CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln','AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu','TGT':'Cys','TGC':'Cys','TGA':'Ter','TGG':'Trp','CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg','AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg','GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'}

# List to hold zero-fold and four-fold sites
site_degen_list = []

# Iterate over each transcript in coordinate file
for transcript in range(0,len(cds_coords)):

    # Get scaffold and id of current transcript
    transcript_id = cds_coords.iloc[transcript]['Name']
    trans_scaff = cds_coords.iloc[transcript]['Scaff']

    # What to do if scaffold is in + direction
    if scaff_orient.loc[scaff_orient.Scaff==trans_scaff, 'Orient'].values[0] == '+':

        # What to do if transcript is in + direction
        if cds_coords.iloc[transcript]['Strand'] == "+":

            # Get sizes of exons from bed12 file into list
            block_sizes = list(cds_coords.iloc[transcript]['BlockSizes'].strip(',').split(','))

            # Get start positions from bed 12 file into list
            block_starts = list(cds_coords.iloc[transcript]['BlockStarts'].strip(',').split(','))

            # Set counters for block number, and position in block
            block = 0
            block_count = 0

            # Get sequence from dictionary and convert to uppercase
            cds_seq = record_dict[transcript_id].upper()
            scaff = cds_coords.iloc[transcript]['Scaff']
            print(transcript_id)

            # Iterate over sequence in intervals of 3 to get codons
            for codon in range(0,len(cds_seq),3):
                
                # set current codon
                current_codon = cds_seq.seq[codon:codon+3]

                # Check each position in codon
                for pos in range(0,len(current_codon)):

                    # Create counter for codon degeneracy
                    # Resets for each position
                    degen = 0

                    # Ignore codon if it contains an 'N'
                    if 'N' in current_codon:
                        degen = 'NA'

                    else:
                        # Replace site with all possible nucleotides
                        for nuc in ['A','C','T','G']:

                            # New codon is set to string of codon up to current position, plus current nucleotide, followed by remainder of the codon
                            new_codon = current_codon[:pos] + nuc + current_codon[pos+1:]
                            # Get original AA for codon
                            og_AA = codon_table[current_codon]
                            # Get AA after replacing nuc and check if the same
                            new_AA = codon_table[new_codon]
                            if og_AA != new_AA:
                                # If not the same, add to counter
                                degen += 1

                    # Convert position from CDS to genome coordinate
                    ## block_starts[block] gives the start position of current exon, with the start of the first exon set to 0
                    ## block_count gives the current position of the nucleotide within the current exon
                    ## cds_coords.iloc[transcript]['Start'] gives the genomic position of the beginning of the gene
                    genome_pos = (int(block_starts[block]) + block_count + cds_coords.iloc[transcript]['Start'])
                    gene_pos = codon + pos + 1
                    # If changing the nucleotide at this position never changed AA, save genomic position and 'Four-fold'
                    if degen == 0:
                        site_degen_list.append([scaff, genome_pos, "Four-fold", current_codon[pos], transcript_id, gene_pos])

                    # If changing to any other nucleotide at this position changed the AA, save genomic position and 'Zero-fold'
                    elif degen == 3:
                        site_degen_list.append([scaff, genome_pos, "Zero-fold", current_codon[pos], transcript_id, gene_pos])

                    # If the current position is the last site in the exon, increase the block number and reset block_count
                    if block_count == int(block_sizes[block]) - 1:
                        block += 1
                        block_count = 0
                    # Otherwise, increase block_count by one
                    else:
                        block_count += 1
        # What do do if scaffold is + and gene strand is -

        elif cds_coords.iloc[transcript]['Strand'] == "-":

            # Get sizes of exons from bed12 file into list and reverse the list order
            block_sizes = list(cds_coords.iloc[transcript]['BlockSizes'].strip(',').split(','))
            block_sizes.reverse()

            # Get start positions from bed 12 file into list and reverse the order
            block_starts = list(cds_coords.iloc[transcript]['BlockStarts'].strip(',').split(','))            
            block_starts.reverse()

            # Set counters for block number, and position in block
            block = 0
            block_count = 0


            # Get CDS seq from dictionary, convert to uppercase 
            cds_seq = record_dict[transcript_id].upper()
            scaff = cds_coords.iloc[transcript]['Scaff']

            print(transcript_id)

            # Iterate over sequence in intervals of 3 to get codons
            for codon in range(0,len(cds_seq),3):

                # set current codon
                current_codon = cds_seq.seq[codon:codon+3]

                # Check each position in codon
                for pos in range(0,len(current_codon)):

                    # Create counter for codon degeneracy
                    degen = 0

                    # Ignore codon if it contains an 'N'
                    if 'N' in current_codon:
                        degen = 'NA'
                    
                    else:
                        # Replace site with all possible nucleotides
                        for nuc in ['A','C','T','G']:

                            # New codon is set to string of codon up to current position, plus current nucleotide, followed by remainder of the codon
                            new_codon = current_codon[:pos] + nuc + current_codon[pos+1:]
                            # Get original AA for codon
                            og_AA = codon_table[current_codon]
                            # Get AA after replacing nuc and check if the same
                            new_AA = codon_table[new_codon]
                            if og_AA != new_AA:
                                # If not the same, add to counter
                                degen += 1

                    # Convert position from CDS to genome coordinate
                    ## block_starts[block] gives the start position of current exon, with the start of the first exon set to 0
                    ## block_count gives the current position of the nucleotide within the current exon
                    ## int(block_sizes[block]) - 1 gives the end position of the exon, because we have reversed the order
                    ## block_count is then subtracted from the end of the exon because we begin at the last position and move towards the begining
                    ## cds_coords.iloc[transcript]['Start'] gives the genomic position of the start of the gene
                    genome_pos = int(block_starts[block]) + (int(block_sizes[block]) - 1) - block_count + cds_coords.iloc[transcript]['Start']
                    gene_pos = codon + pos + 1
                    # If changing the nucleotide at this position never changed AA, save genomic position and 'Four-fold'
                    if degen == 0:
                        site_degen_list.append([scaff, genome_pos, "Four-fold", Seq(current_codon[pos]).reverse_complement(), transcript_id, gene_pos])

                    # If changing to any other nucleotide at this position changed the AA, save genomic position and 'Zero-fold'
                    elif degen == 3:
                        site_degen_list.append([scaff, genome_pos, "Zero-fold", Seq(current_codon[pos]).reverse_complement(), transcript_id, gene_pos])

                    # If the current position is the last site in the exon, increase the block number and reset block_count
                    if block_count == int(block_sizes[block]) - 1:
                        block += 1
                        block_count = 0
                    # Otherwise, increase block_count by one
                    else:
                        block_count += 1

    # What to do if the scaffold orientation is - 
    # In this case, the actions are reversed for genes with + strand vs. - strand
    elif scaff_orient.loc[scaff_orient.Scaff==trans_scaff, 'Orient'].values[0] == '-':

        # What to do if gene strand is -
        if cds_coords.iloc[transcript]['Strand'] == "-":

            # Get sizes of exons from bed12 file into list
            block_sizes = list(cds_coords.iloc[transcript]['BlockSizes'].strip(',').split(','))

            # Get start positions from bed 12 file into list
            block_starts = list(cds_coords.iloc[transcript]['BlockStarts'].strip(',').split(','))

            # Set counters for block number, and position in block
            block = 0
            block_count = 0



            # Get sequence from dictionary and convert to uppercase
            cds_seq = record_dict[transcript_id].upper()
            scaff = cds_coords.iloc[transcript]['Scaff']

            print(transcript_id)

            # Iterate over sequence in intervals of 3 to get codons
            for codon in range(0,len(cds_seq),3):

                # set current codon
                current_codon = cds_seq.seq[codon:codon+3]

                # Check each position in codon
                for pos in range(0,len(current_codon)):

                    # Create counter for codon degeneracy
                    degen = 0

                    # Ignore codon if it contains an 'N'
                    if 'N' in current_codon:
                        degen = 'NA'

                    else:
                        # Replace site with all possible nucleotides
                        for nuc in ['A','C','T','G']:

                            # New codon is set to string of codon up to current position, plus current nucleotide, followed by remainder of the codon
                            new_codon = current_codon[:pos] + nuc + current_codon[pos+1:]
                            # Get original AA for codon
                            og_AA = codon_table[current_codon]
                            # Get AA after replacing nuc and check if the same
                            new_AA = codon_table[new_codon]
                            if og_AA != new_AA:
                                # If not the same, add to counter
                                degen += 1

                    # Convert position from CDS to genome coordinate
                    ## block_starts[block] gives the start position of current exon, with the start of the first exon set to 0
                    ## block_count gives the current position of the nucleotide within the current exon
                    ## cds_coords.iloc[transcript]['Start'] gives the genomic position of the beginning of the gene
                    genome_pos = (int(block_starts[block]) + block_count + cds_coords.iloc[transcript]['Start'])
                    gene_pos = codon + pos + 1
                    # If changing the nucleotide at this position never changed AA, save genomic position and 'Four-fold'
                    if degen == 0:
                        site_degen_list.append([scaff, genome_pos, "Four-fold", current_codon[pos], transcript_id, gene_pos])
                    # If changing to any other nucleotide at this position changed the AA, save genomic position and 'Zero-fold'
                    elif degen == 3:
                        site_degen_list.append([scaff, genome_pos, "Zero-fold", current_codon[pos], transcript_id, gene_pos])

                    # If the current position is the last site in the exon, increase the block number and reset block_count
                    if block_count == int(block_sizes[block]) - 1:
                        block += 1
                        block_count = 0
                    # Otherwise, increase block_count by one
                    else:
                        block_count += 1
                    
        # What to do if scaffold orientation is - and gene strand is +
        elif cds_coords.iloc[transcript]['Strand'] == "+":

            # Get sizes of exons from bed12 file into list and reverse the list order
            block_sizes = list(cds_coords.iloc[transcript]['BlockSizes'].strip(',').split(','))
            block_sizes.reverse()

            # Get start positions from bed 12 file into list and reverse the order
            block_starts = list(cds_coords.iloc[transcript]['BlockStarts'].strip(',').split(','))
            block_starts.reverse()

            # Set counters for block number, and position in block
            block = 0
            block_count = 0

            # Get CDS seq from dictionary, convert to uppercase 
            cds_seq = record_dict[transcript_id].upper()
            scaff = cds_coords.iloc[transcript]['Scaff']
            print(transcript_id)

            # Iterate over sequence in intervals of 3 to get codons
            for codon in range(0,len(cds_seq),3):

                # set current codon
                current_codon = cds_seq.seq[codon:codon+3]

                # Check each position in codon
                for pos in range(0,len(current_codon)):

                    # Create counter for codon degeneracy
                    degen = 0

                    # Ignore codon if it contains an 'N'
                    if 'N' in current_codon:
                        degen = 'NA'
                    else:

                        # Replace site with all possible nucleotides
                        for nuc in ['A','C','T','G']:

                            # New codon is set to string of codon up to current position, plus current nucleotide, followed by remainder of the codon
                            new_codon = current_codon[:pos] + nuc + current_codon[pos+1:]

                            # Get original AA for codon
                            og_AA = codon_table[current_codon]

                            # Get AA after replacing nuc and check if the same
                            new_AA = codon_table[new_codon]
                            if og_AA != new_AA:
                                # If not the same, add to counter
                                degen += 1

                    # Convert position from CDS to genome coordinate
                    ## block_starts[block] gives the start position of current exon, with the start of the first exon set to 0
                    ## block_count gives the current position of the nucleotide within the current exon
                    ## int(block_sizes[block]) - 1 gives the end position of the exon, because we have reversed the order
                    ## block_count is then subtracted from the end of the exon because we begin at the last position and move towards the begining
                    ## cds_coords.iloc[transcript]['Start'] gives the genomic position of the start of the gene
                    genome_pos = int(block_starts[block]) + (int(block_sizes[block]) - 1) - block_count + cds_coords.iloc[transcript]['Start']
                    gene_pos = codon + pos + 1
                    # If changing the nucleotide at this position never changed AA, save genomic position and 'Four-fold'
                    if degen == 0:
                        site_degen_list.append([scaff, genome_pos, "Four-fold", Seq(current_codon[pos]).reverse_complement(), transcript_id, gene_pos])

                    # If changing to any other nucleotide at this position changed the AA, save genomic position and 'Zero-fold'
                    elif degen == 3:
                        site_degen_list.append([scaff, genome_pos, "Zero-fold", Seq(current_codon[pos]).reverse_complement(), transcript_id, gene_pos])

                    # If the current position is the last site in the exon, increase the block number and reset block_count
                    if block_count == int(block_sizes[block]) - 1:
                        block += 1
                        block_count = 0

                    # Otherwise, increase block_count by one
                    else:
                        block_count += 1

# Create new dataframe from site_degen_list
out_df = pd.DataFrame(site_degen_list, columns=['Scaff','Pos','Degen','Anc_nuc', 'Trans_id', 'Gene_pos'])


# Add 1 to all positions to end up in proper genome position
out_df['Pos'] = out_df['Pos'] + 1

# Positions are out of order because of working with reverse transcripts
# Sort values by scaffold and position
out_df.sort_values(by = ['Scaff','Pos'], inplace = True, ignore_index=True)

# Write dataframe to file
out_df.to_csv(out_file, sep = '\t', index = False)

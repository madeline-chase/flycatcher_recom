import pandas as pd

files = pd.read_table('aln_gene_names.txt', names = ['file_name'])
localrules: all

rule all:
    input:
        aligned = expand('results/dn_ds/aln/fAlb15/prank/{gene}.clustal_tree.aln.best.fas', gene = files['file_name'])

rule clustal:
    input: 
        seq = 'intermediate/dn_ds/concat_seqs/fAlb15/{gene}.fasta'
    output:
        clustal_aln = 'results/dn_ds/aln/fAlb15/clustal/{gene}.clustal_aln_out.fas',
        clustal_tree = 'intermediate/dn_ds/concat_seqs/fAlb15/{gene}.dnd'
    shell:
        """
        clustalw2 -INFILE={input.seq} -ALIGN -OUTFILE={output.clustal_aln} -OUTPUT=FASTA
        """

rule prank:
    input:
        seq = 'intermediate/dn_ds/concat_seqs/fAlb15/{gene}.fasta',
        tree = 'intermediate/dn_ds/concat_seqs/fAlb15/{gene}.dnd'
    params:
        prefix= 'results/dn_ds/aln/fAlb15/prank/{gene}.clustal_tree.aln'
    output:
        prank_aln = 'results/dn_ds/aln/fAlb15/prank/{gene}.clustal_tree.aln.best.fas'
    shell:
        """
        prank -d={input.seq} -o={params.prefix} -F -codon -t={input.tree} -once
        """
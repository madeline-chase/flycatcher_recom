## Snakemake workflow to run recombination rate estimation workflow using LDhelmet
configfile: "config.yaml"
localrules: all

import pandas as pd

scaffolds = pd.read_table(config['scaff_file'], names = ['chr_name','scaff_name', 'category'])

rule all:
    input:
        sf_out = expand('results/sf_out/{spec}/no_cpg_filt/{spec}.{scaff}.sf2.ug.anc_fix_rmvd.no_cpg_filt.out.txt', spec = config['spec'], scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'])

## Get count data for polarized sites
rule get_count_data:
    input:
        vcf = '../allsites_vc/data/vcf/non_z_scaff/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.vcf.gz',
        polar_sites = '../allsites_vc/data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.non_z_scaff.var_only_counts.groups_combined.txt',
        samps = '../allsites_vc/data/samples/{spec}.all.txt'
    output:
        counts = 'data/counts/no_cpg_filt/{spec}.{scaff}.polar_sites.no_cpg_filt.counts.frq.count'
    params:
        prefix = 'data/counts/no_cpg_filt/{spec}.{scaff}.polar_sites.no_cpg_filt.counts'
    shell:
        """
        vcftools --gzvcf {input.vcf} --counts --positions {input.polar_sites} --keep {input.samps} --out {params.prefix}
        """

## Split count data output
rule split_count_data:
    input:
        counts = 'data/counts/no_cpg_filt/{spec}.{scaff}.polar_sites.no_cpg_filt.counts.frq.count'
    output:
        counts_split = 'data/counts/no_cpg_filt/{spec}.{scaff}.polar_sites.no_cpg_filt.counts.split.txt'
    shell:
        """
        awk -v OFS='\t' '{{if(NR>1) {{split($5,a,":") split($6,b,":"); print $1,$2,$3,$4,a[1],a[2],b[1],b[2]}}}}' {input.counts} > {output.counts_split}
        """

## Reformat to get counts of derived alleles for SF2 input
## Only output sites that are not fixed for ancestral allele
rule get_sf_data:
    input:
        polar_sites = '../allsites_vc/data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.non_z_scaff.var_only_counts.groups_combined.txt',
        split_counts = 'data/counts/no_cpg_filt/{spec}.{scaff}.polar_sites.no_cpg_filt.counts.split.txt',
        sf_header = 'data/sf_header.txt'
    output:
        sf_count_input = 'data/counts/sf_in/no_cpg_filt/{spec}.{scaff}.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.txt'
    shell:
        """
        paste {input.polar_sites} {input.split_counts} | awk -v OFS='\t' '{{if($3==$8) print $5,$11,$7,"0"; else if($3==$10) print $5,$9,$7,"0"}}' | awk '{{if($2>0) print}}' | cat {input.sf_header} -  > {output.sf_count_input}
        """

## Get grid file for SF2 input
## If no sites present for scaffold just create empty file
rule get_grid_file:
    input:
        sf_count_input = 'data/counts/sf_in/no_cpg_filt/{spec}.{scaff}.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.txt'
    output:
        sf_grid_file = 'data/counts/sf_in/no_cpg_filt/{spec}.{scaff}.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.grid_file.txt'
    run:
        if len(pd.read_csv(input.sf_count_input)) > 0:
            shell("grep -v 'folded' {input.sf_count_input} | awk '{{print $1}}' > {output.sf_grid_file}")
        else:
            shell("touch {output.sf_grid_file}")

## Combine all scaffolds for whole genome input for background SFS
rule get_whole_genome_input:
    input:
        counts = expand('data/counts/sf_in/no_cpg_filt/{{spec}}.{scaff}.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.txt', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        sf_header = 'data/sf_header.txt'
    output:
        whole_gen_counts = 'data/counts/sf_in/no_cpg_filt/{spec}.whole_genome.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.txt'
    shell:
        """
        cat {input.counts} | grep -v 'folded' | cat {input.sf_header} - > {output.whole_gen_counts}
        """

## Get whole genome SFS
rule get_whole_genome_sf:
    input:
        whole_gen_counts = 'data/counts/sf_in/no_cpg_filt/{spec}.whole_genome.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.txt'
    output:
        whole_gen_spec = 'data/sf_spec/no_cpg_filt/{spec}.whole_genome.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.sf_spec'
    shell:
        """
        scripts/SF2/SweepFinder2 -f {input.whole_gen_counts} {output.whole_gen_spec}
        """

## Run SF2 for all scaffolds, if enough sites present on scaffold (>15)
rule run_sf:
    input:
        sf_counts = 'data/counts/sf_in/no_cpg_filt/{spec}.{scaff}.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.txt',
        whole_gen_spec = 'data/sf_spec/no_cpg_filt/{spec}.whole_genome.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.sf_spec',
        sf_grid_file = 'data/counts/sf_in/no_cpg_filt/{spec}.{scaff}.polar_sites.sf2_format.anc_fix_rmvd.no_cpg_filt.grid_file.txt'
    output:
        sf_out = 'results/sf_out/{spec}/no_cpg_filt/{spec}.{scaff}.sf2.ug.anc_fix_rmvd.no_cpg_filt.out.txt'
    run:
        if len(pd.read_csv(input.sf_counts)) > 15:
            shell("scripts/SF2/SweepFinder2 -lu {input.sf_grid_file} {input.sf_counts} {input.whole_gen_spec} {output.sf_out}")
        else:
            shell("touch {output.sf_out}")


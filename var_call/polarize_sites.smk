configfile: 'polarize_snps_config.yaml'
localrules: all

import pandas as pd

scaffolds = pd.read_table(config['scaff_file'], names = ['chr_name','scaff_name', 'category'])

rule all:
    input:
        anc_fasta = 'data/ref/anc_fasta/fAlb15_MTmask_mtfZan.chrom_scaffs_only.noncall_masked.snp_masked.anc_allele.no_cpg_filt.fasta'


## Get counts for all sites in hyp, removing missing sites
rule hyp_counts:
    input:
        vcf = 'data/vcf/hyp/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only.vcf.gz',
        missing = 'data/vcf/{scaff_type}/missing_sites/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.all_spec_comb.missing_sites_10_perc.txt'
    output:
        counts = 'data/counts/missing_filt/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.frq.count'
    params:
        prefix = 'data/counts/missing_filt/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}'
    shell:
        """
        vcftools --gzvcf {input.vcf} --exclude-positions {input.missing} --out {params.prefix} --counts
        """

## Subset counts for variant sites only for hyp
rule get_var_hyp:
    input:
        counts = 'data/counts/missing_filt/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.frq.count'
    output:
        var_counts = 'data/counts/var_only/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.txt'
    shell:
        """
        awk '{{if($3==2) print}}' {input.counts} > {output.var_counts}
        """

## Split count file for hyp
rule split_counts_hyp:
    input:
        counts = 'data/counts/var_only/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.txt'
    output:
        split_counts =  'data/counts/var_only/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.split.txt'
    shell:
        """
        awk -v OFS='\t' '{{split($5,a,":") split($6,b,":"); print $1,$2,$3,$4,a[1],a[2],b[1],b[2]}}' {input.counts} > {output.split_counts}
        """

## Get counts for all sites in collared+pied and red-breasted+taig, removing missing sites
rule af_counts:
    input:
        vcf = 'data/vcf/{scaff_type}/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.vcf.gz',
        missing = 'data/vcf/{scaff_type}/missing_sites/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.all_spec_comb.missing_sites_10_perc.txt',
        samp = 'data/samples/{group_name}.all.txt'
    output:
        counts = 'data/counts/missing_filt/all_sites.{scaff}.all_filt.missing_rmvd.{group_name}.{scaff_type}.frq.count'
    params:
        prefix = 'data/counts/missing_filt/all_sites.{scaff}.all_filt.missing_rmvd.{group_name}.{scaff_type}'
    shell:
        """
        vcftools --gzvcf {input.vcf} --exclude-positions {input.missing} --out {params.prefix} --counts --keep {input.samp}
        """

## Subset counts for only variable sites
rule get_var_sites:
    input:
        counts = 'data/counts/missing_filt/all_sites.{scaff}.all_filt.missing_rmvd.{group_name}.{scaff_type}.frq.count'
    output:
        var_counts = 'data/counts/var_only/all_sites.{scaff}.all_filt.missing_rmvd.{group_name}.{scaff_type}.var_only_counts.txt'
    shell:
        """
        awk '{{if($3==2) print}}' {input.counts} > {output.var_counts}
        """

## Split count files
rule split_counts:
    input:
        counts = 'data/counts/var_only/all_sites.{scaff}.all_filt.missing_rmvd.{group_name}.{scaff_type}.var_only_counts.txt'
    output:
        split_counts  = 'data/counts/var_only/all_sites.{scaff}.all_filt.missing_rmvd.{group_name}.{scaff_type}.var_only_counts.split.txt'
    shell:
        """
        awk -v OFS='\t' '{{split($5,a,":") split($6,b,":"); print $1,$2,$3,$4,a[1],a[2],b[1],b[2]}}' {input.counts} > {output.split_counts}
        """

## Write bed format file for snp positions
rule snps_to_bed:
    input:
        split_counts = 'data/counts/var_only/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.split.txt'
    output:
        bed_counts = 'data/counts/var_only/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.split.bed'
    shell:
        """
        awk -v OFS='\t' '{{print $1,$2-1,$2}}' {input.split_counts} > {output.bed_counts}
        """

## Combine all SNPs for each scaffold
rule combine_snps_auto:
    input:
        bed_counts = expand('data/counts/var_only/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.split.bed', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'], scaff_type = 'non_z_scaff')
    output:
        snps_merged = 'data/counts/var_only/all_sites.hyp_only.all_auto_scaff.all_filt.missing_rmvd.var_only_counts.split.bed'
    shell:
        """
        cat {input.bed_counts} | sort -k1,1 -k2n,2 > {output.snps_merged}
        """

## Combine var site count data from three groups into one file
rule combine_three_group_counts:
    input:
        hyp_counts = 'data/counts/var_only/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.split.txt',
        coll_pied_counts = 'data/counts/var_only/all_sites.{scaff}.all_filt.missing_rmvd.coll_pied.{scaff_type}.var_only_counts.split.txt',
        par_taig_counts = 'data/counts/var_only/all_sites.{scaff}.all_filt.missing_rmvd.par_taig.{scaff_type}.var_only_counts.split.txt'
    output:
        counts_combined = 'data/counts/groups_combined/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.groups_combined.txt'
    shell:
        """
        paste {input.hyp_counts} {input.coll_pied_counts} {input.par_taig_counts} > {output.counts_combined}
        """

## Identify sites where 2 of 3 groups fixed for same allele
rule identify_ancestral_allele:
    input:
        counts_combined = 'data/counts/groups_combined/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.groups_combined.txt'
    output:
        anc_allele = 'data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.groups_combined.txt'
    shell:
        """
        python scripts/get_anc_allele.py {input.counts_combined} {output.anc_allele}
        """

## Get bed formatted file for ancestral sites
rule anc_allele_bed:
    input:
        anc_allele = 'data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.groups_combined.txt'
    output:
        anc_bed = 'data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.groups_combined.bed'
    shell:
        """
        awk -v OFS='\t' '{{print $1,$2-1,$2,$3}}' {input.anc_allele} > {output.anc_bed}
        """

## Combine ancestral sites for all scaffolds
rule combine_anc_auto:
    input:
        anc_alleles = expand('data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.var_only_counts.groups_combined.bed',scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'], scaff_type = 'non_z_scaff')
    output:
        anc_merged = 'data/counts/polar_sites/no_cpg_filt/polar_sites.all_auto_scaff.all_filt.missing_rmvd.var_only_counts.groups_combined.bed'
    shell:
        """
        cat {input.anc_alleles} | sort -k1,1 -k2n,2 > {output.anc_merged}
        """

## Using counts, convert sites to bed format and merge into contiguous regions
## these then give the regions where we have callable data
rule counts_to_bed:
    input:
        counts = 'data/counts/missing_filt/all_sites.hyp_only.{scaff}.all_filt.missing_rmvd.{scaff_type}.frq.count'
    output:
        bed_counts = 'data/callable_data/called_regions/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.callable_sites.bed'
    shell:
        """
        awk -v OFS='\t' '{{if(NR>1) {{print $1,$2-1,$2}}}}' {input.counts} > {output.bed_counts}
        """

## Merge counts to get called regions
rule get_called_regions:
    input:
        bed_counts = 'data/callable_data/called_regions/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.callable_sites.bed'
    output:
        bed_merged = 'data/callable_data/called_regions/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.callable_sites.merged.bed'
    shell:
        """
        bedtools merge -i {input.bed_counts} > {output.bed_merged}
        """

## Use bedtools to subtract the callabe sites from the genome file
rule get_noncall_regions:
    input:
        bed_merged = 'data/callable_data/called_regions/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.callable_sites.merged.bed',
        genome_file = 'data/ref/fAlb15_MTmas_mtfZan.genomeFile.chrom.all.20140121.bed'
    output:
        noncall_bed = 'data/callable_data/noncall_regions/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.noncall_sites.bed'
    shell:
        """
        grep {wildcards.scaff} {input.genome_file} | bedtools subtract -a stdin -b {input.bed_merged} > {output.noncall_bed}
        """

## Combine noncall regions for all scaffolds
rule combine_noncall_regions_auto:
    input:
        noncall_beds = expand('data/callable_data/noncall_regions/all_sites.{scaff}.all_filt.missing_rmvd.{scaff_type}.noncall_sites.bed', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'], scaff_type = 'non_z_scaff')
    output:
        noncall_merged = 'data/callable_data/noncall_regions/all_sites.all_auto_scaff.all_filt.missing_rmvd.noncall_sites.bed'
    shell:
        """
        cat {input.noncall_beds} | sort -k1,1 -k2n,2 > {output.noncall_merged}
        """

## Use bedtools maskfasta to mask noncallable sites from reference genome
rule mask_noncall_sites:
    input:
        noncall_bed ='data/callable_data/noncall_regions/all_sites.all_auto_scaff.all_filt.missing_rmvd.noncall_sites.bed',
        fasta = 'data/ref/fAlb15_MTmask_mtfZan.chrom_scaffs_only.fasta'
    output:
        masked_fasta = 'data/ref/fAlb15_MTmask_mtfZan.chrom_scaffs_only.noncall_masked.fasta'
    shell:
        """
        bedtools maskfasta -fi {input.fasta} -bed {input.noncall_bed} -fo {output.masked_fasta}
        """

## Use maskfasta to mask all variable sites in ref genome
rule mask_snps:
    input:
        snp_data = 'data/counts/var_only/all_sites.hyp_only.all_auto_scaff.all_filt.missing_rmvd.var_only_counts.split.bed',
        masked_fasta = 'data/ref/fAlb15_MTmask_mtfZan.chrom_scaffs_only.noncall_masked.fasta'
    output:
        snp_masked_fasta = 'data/ref/fAlb15_MTmask_mtfZan.chrom_scaffs_only.noncall_masked.snp_masked.fasta'
    shell:
        """
        bedtools maskfasta -fi {input.masked_fasta} -bed {input.snp_data} -fo {output.snp_masked_fasta}
        """

## Use script to replace masked sites with ancestral allele where identified
rule replace_ancestral_allele:
    input:
        anc_merged = 'data/counts/polar_sites/no_cpg_filt/polar_sites.all_auto_scaff.all_filt.missing_rmvd.var_only_counts.groups_combined.bed',
        snp_masked_fasta = 'data/ref/fAlb15_MTmask_mtfZan.chrom_scaffs_only.noncall_masked.snp_masked.fasta'
    output:
        anc_fasta = 'data/ref/anc_fasta/fAlb15_MTmask_mtfZan.chrom_scaffs_only.noncall_masked.snp_masked.anc_allele.no_cpg_filt.fasta'
    shell:
        """
        cp {input.snp_masked_fasta} {output.anc_fasta}
        python scripts/replace_fasta_bases.py {input.anc_merged} {output.anc_fasta} 
        """


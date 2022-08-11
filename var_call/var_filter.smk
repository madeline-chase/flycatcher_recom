configfile: 'filter_vars_config.yaml'
localrules: all, combine_z_stats, combine_non_z_stats

import pandas as pd

scaffolds = pd.read_table(config['scaff_file'], names = ['chr_name','scaff_name', 'category'])

rule all:
    input:
        non_z_missings_stats = 'data/stats/non_z_scaff/missing_filt/all_stats_file_names.txt',
        non_z_CM_stats='data/stats/non_z_scaff/collapsed_filt/all_stats_file_names.txt'

## Get stats for unfiltered vcfs
rule bcf_stats:
    input:
        vcf = 'data/vcf/all_sites.{scaff}.vcf.gz'
    output:
        stats = 'data/stats/unfiltered/all_sites.{scaff}.vcf_stats.txt'
    shell:
        """
        bcftools stats {input.vcf} > {output.stats}
        """

## Remove indels and biallelic sites
rule remove_indels:
    input:
        vcf = 'data/vcf/unfilt/all_sites.{scaff}.vcf.gz',
        stats = 'data/stats/unfiltered/all_sites.{scaff}.vcf_stats.txt'
    output:
        vcf_out = 'data/vcf/indel_rmvd/all_sites.{scaff}.indel_rmvd.max_allele2.vcf.gz'
    shell:
        """
        vcftools --gzvcf {input.vcf} --remove-indels --recode --max-alleles 2 --recode-INFO-all --stdout | gzip > {output.vcf_out}
        """

## Apply hard filter cutoffs to vcfs
rule gatk_hard_filt:
    input: 
        vcf = 'data/vcf/bgzip/all_sites.{scaff}.indel_rmvd.max_allele2.vcf.gz',
        ref = '/proj/sllstore2017033/nobackup/b2017012_nobackup/reference/fAlb15_MTmask_mtfZan.fasta'
    output: 
        vcf_out = 'data/vcf/hard_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.vcf.gz'
    shell:
        """
        module load bioinfo-tools GATK/4.1.4.1
        gatk VariantFiltration -R {input.ref} -V {input.vcf} -O {output.vcf_out} --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || StrandOddsRatio > 3 || ReadPosRankSum < -8.0" --filter-name "hard_filt" -L {wildcards.scaff}
        """

## Remove sites that don't pass hard filters
rule remove_filt:
    input:
        vcf = 'data/vcf/hard_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.vcf.gz'
    output:
        vcf_out = 'data/vcf/hard_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.vcf.gz'
    shell:
        """
        vcftools --gzvcf {input.vcf} --remove-filtered-all --recode --recode-INFO-all --stdout | gzip > {output.vcf_out}
        """

## Remove sites overlapping with repeats
rule remove_repeats:
    input:
        vcf = 'data/vcf/hard_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.vcf.gz',
        repeats = 'data/repeats/fAlb15.fa.out.bed'
    output:
        vcf_out = 'data/vcf/repeat_masked/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.vcf.gz'
    shell:
        """
        vcftools --gzvcf {input.vcf} --exclude-bed {input.repeats} --recode --recode-INFO-all --stdout | gzip > {output.vcf_out}   
        """

## Apply GT thresholds to hyperythra sample
rule hyp_gt_filt:
    input:
        vcf = 'data/vcf/repeat_masked/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.vcf.gz'
    output:
        vcf_out = 'data/vcf/hyp/gt_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.minDP.maxDP.hyp_only.vcf.gz'
    shell:
        """
        vcftools --gzvcf {input.vcf} --minDP 5 --maxDP 200 --recode --recode-INFO-all --indv F_hyperythra_X --stdout | gzip > {output.vcf_out}
        """

## Remove collapsed regions from hyp sample
rule collapsed_filt_hyp:
    input:
        vcf = 'data/vcf/hyp/gt_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.minDP.maxDP.hyp_only.vcf.gz',
        collapsed = 'data/vcf/collapsed_reg/merged_bed_species/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.c7.e3.n5.d1000.x2000.all_merged.m5000.bed'
    output:
        vcf_out = 'data/vcf/hyp/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only.vcf.gz'
    shell:
        """
        vcftools --gzvcf {input.vcf} --exclude-bed {input.collapsed} --recode --recode-INFO-all --stdout | gzip > {output.vcf_out}
        """

## Get allele counts for hyp    
rule hyp_af_count:
    input:
        vcf = 'data/vcf/hyp/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only.vcf.gz'
    output:
        counts = 'data/vcf/hyp/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only.frq.count'
    params:
        prefix = 'data/vcf/hyp/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only'
    shell:
        """
        vcftools --gzvcf {input.vcf} --counts --out {params.prefix}
        """

## Get missing sites for hyp from allele counts
rule hyp_missing:
    input:
        counts = 'data/vcf/hyp/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only.frq.count'
    output:
        missing = 'data/vcf/hyp/missing_sites/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only.missing_sites.txt'
    shell:
        """
        awk -v OFS='\t' '{{if(NR>1) {{if($4==0) print $1,$2}}}}' {input.counts} > {output.missing}
        """

## Apply GT thresholds to all other species
rule non_z_gt_filt:
    input: 
        vcf = 'data/vcf/repeat_masked/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.vcf.gz'
    output:
        vcf_out = 'data/vcf/non_z_scaff/gt_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.minGQ.minDP.maxDP.vcf.gz'
    shell:
        """
        vcftools --gzvcf {input.vcf} --minGQ 30 --minDP 5 --maxDP 200 --recode --recode-INFO-all --stdout | gzip > {output.vcf_out}
        """

## Get vcf stats after filter
rule post_filt_stats_non_z:
    input:
        vcf = 'data/vcf/non_z_scaff/gt_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.minGQ.minDP.maxDP.vcf.gz',
    output:
        stats = 'data/stats/non_z_scaff/gt_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.minGQ.minDP.maxDP.stats.txt'
    shell:
        """
        bcftools stats {input.vcf} > {output.stats}
        """

## Remove sites overlapping with collapsed regions
rule filter_collapsed_non_z:
    input:
        vcf = 'data/vcf/non_z_scaff/gt_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.minGQ.minDP.maxDP.vcf.gz',
        stats = 'data/stats/non_z_scaff/gt_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.minGQ.minDP.maxDP.stats.txt',
        collapsed = 'data/vcf/collapsed_reg/merged_bed_species/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.c7.e3.n5.d1000.x2000.all_merged.m5000.bed'
    output:
        vcf_out = 'data/vcf/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.vcf.gz'
    shell:
        """
        vcftools --gzvcf {input.vcf} --exclude-bed {input.collapsed} --recode --recode-INFO-all --stdout | gzip > {output.vcf_out}
        """

## Get allele counts for each species separately 
rule get_af_by_species_non_z:
    input:
        vcf = 'data/vcf/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.vcf.gz',
        samples = 'data/samples/{spec}.all.txt'
    output:
        counts = 'data/vcf/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.{spec}.frq.count'
    params:
        prefix = 'data/vcf/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.{spec}'
    shell:
        """
        vcftools --gzvcf {input.vcf} --counts --keep {input.samples} --out {params.prefix}
        """

## Get missing sites for each species with at least 10% missing data
rule find_missing_sites_non_z:
    input:
        counts = 'data/vcf/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.{spec}.frq.count',
        samples = 'data/samples/{spec}.all.txt'
    output:
        missing_sites = 'data/vcf/non_z_scaff/missing_sites/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.{spec}.missing_sites_10_perc.txt'
    shell:
        """
        ind=`cat {input.samples} | wc -l`
        awk -v OFS='\t' -v var=$ind '{{if(NR>1) {{if($4<(2*var)*0.9) print $1,$2}}}}' {input.counts} > {output.missing_sites}

        """

## Combine missing sites for all species
rule combine_missing_sites_non_z:
    input:
        missing_sites = expand('data/vcf/non_z_scaff/missing_sites/all_sites.{{scaff}}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.{spec}.missing_sites_10_perc.txt', spec = config['spec']),
        hyp_missing = 'data/vcf/hyp/missing_sites/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minDP.maxDP.hyp_only.missing_sites.txt'
    output:
        missing_combined = 'data/vcf/non_z_scaff/missing_sites/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.all_spec_comb.missing_sites_10_perc.txt'
    shell:
        """
        cat {input.missing_sites} {input.hyp_missing} | sort -k1,1 -k2n,2 | uniq > {output.missing_combined}
        """

## Filter missing sites
rule remove_missing_sites_non_z:
    input:
        vcf = 'data/vcf/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.vcf.gz',
        missing_sites = 'data/vcf/non_z_scaff/missing_sites/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.all_spec_comb.missing_sites_10_perc.txt'
    output:
        vcf_out = 'data/vcf/non_z_scaff/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.vcf.gz',
    shell:
        """
        vcftools --gzvcf {input.vcf} --exclude-positions {input.missing_sites} --recode --recode-INFO-all --stdout | gzip > {output.vcf_out}
        """

## Remove collapsed regions
rule collapsed_filt_stats_non_z:
    input:
        vcf = 'data/vcf/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.vcf.gz'
    output:
        stats = 'data/stats/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.stats.txt'
    shell:
        """
        bcftools stats {input.vcf} > {output.stats}
        """

## Get stats for final filtered vcfs
rule missing_get_filt_stats_non_z:
    input:
        vcf = 'data/vcf/non_z_scaff/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.vcf.gz'
    output:
        stats = 'data/stats/non_z_scaff/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.stats.txt'
    shell:
        """
        bcftools stats {input.vcf} > {output.stats}
        """

## Combine stats files for all scaffolds
rule combine_non_z_stats:
    input:
        missing_stats = expand('data/stats/non_z_scaff/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.stats.txt', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        CM_stats = expand('data/stats/non_z_scaff/collapsed_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.stats.txt', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'])
    output:
        missing_stats_file = 'data/stats/non_z_scaff/missing_filt/all_stats_file_names.txt',
        CM_stats_file = 'data/stats/non_z_scaff/collapsed_filt/all_stats_file_names.txt'
    shell:
        """
        echo {input.missing_stats} > {output.missing_stats_file}
        echo {input.CM_stats} > {output.CM_stats_file}
        """


   

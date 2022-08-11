## Snakemake workflow to run recombination rate estimation workflow using LDhelmet
configfile: "config.yaml"
localrules: all,post_to_text

import pandas as pd

scaffolds = pd.read_table(config['scaff_file'], names = ['chr_name','scaff_name', 'category'])


rule all:
    input:
        ldhelm_out = expand( "intermediate/ldhelm/rho/{spec}/no_cpg_filt/{spec}.{scaff}.{phase}.{iter}.rjmcmc.txt", scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'], spec = config["spec"], phase = config["phase_attempts"], iter = config["rjmcmc_its"])

## Convert vcf to plink by species and remove singletons
rule vcf_to_plink:
    input: 
        vcf = '../allsites_vc/data/vcf/non_z_scaff/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.vcf.gz',
        samples = '../allsites_vc/data/samples/{spec}.all.txt'
    output:
        ped = "intermediate/ldhelm/shapeit_in/{spec}/cpg_filt/{spec}.{scaff}.ped",
        map = "intermediate/ldhelm/shapeit_in/{spec}/cpg_filt/{spec}.{scaff}.map"
    params:
        vcftools_pre = "intermediate/ldhelm/shapeit_in/{spec}/{spec}.{scaff}"
    shell:
        """
        vcftools --gzvcf {input.vcf} --plink --out {params.vcftools_pre} --keep {input.samples} --mac 2
        """

## Print anc priors for each scaff from polarized sites + total sites
rule get_anc_priors:
    input:
        anc_alleles = '../allsites_vc/data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.non_z_scaff.var_only_counts.groups_combined.txt',
        map = "intermediate/ldhelm/shapeit_in/{spec}/cpg_filt/{spec}.{scaff}.map"
    output:
        anc_priors = 'intermediate/ldhelm/anc_priors/{spec}/{spec}.{scaff}.anc_priors.no_cpg_filt.zero_pos.txt'
    shell:
        """
        python scripts/print_anc_priors.py {input.map} {input.anc_alleles} > {output.anc_priors}
        """

## Run shapeit for two iterations
rule shapeit:
    input:
        ped = "intermediate/ldhelm/shapeit_in/{spec}/cpg_filt/{spec}.{scaff}.ped",
        map = "intermediate/ldhelm/shapeit_in/{spec}/cpg_filt/{spec}.{scaff}.map"
    output:
        hap = "intermediate/ldhelm/shapeit_out/{spec}/{spec}.{scaff}.{phase}.haps",
        samp = "intermediate/ldhelm/shapeit_out/{spec}/{spec}.{scaff}.{phase}.sample"
    params:
        shapeit_prefix = "intermediate/ldhelm/shapeit_out/{spec}/{spec}.{scaff}.{phase}",
        pop_size = config['pop_size'],
        rho = config['rho'],
        states = config['states'],
        window = config['window']
    shell:
        """
        module load bioinfo-tools SHAPEIT/v2.r904
        shapeit --input-ped {input.ped} {input.map} --window {params.window} --effective-size {params.pop_size} --rho {params.rho} -O {params.shapeit_prefix} --seed $RANDOM --force --states {params.states}
        """

## Convert output from shapeit into ldhelmet input format and subsample haps to 25 specified individuals
rule shapeit2ldhelm:
    input:
        hap = "intermediate/ldhelm/shapeit_out/{spec}/cpg_filt/{spec}.{scaff}.{phase}.haps",
        samp = "intermediate/ldhelm/shapeit_out/{spec}/cpg_filt/{spec}.{scaff}.{phase}.sample",
        keep = "data/{spec}.top_25_samples.txt"
    output:
        ldhelm_seq = "intermediate/ldhelm/ldhelm_in/{spec}/cpg_filt/{spec}.{scaff}.{phase}.seq",
        pos = "intermediate/ldhelm/ldhelm_in/{spec}/cpg_filt/{spec}.{scaff}.{phase}.zero_base.pos",
    shell:
        """
        python scripts/shapeit2ldhelm.py --hap {input.hap} --samp {input.samp} --keep {input.keep} > {output.ldhelm_seq}
        awk '{{print $3-1}}' {input.hap} > {output.pos}
        """

## run ldhelm find_confs for all scaffs combined
rule ldhelm_conf:
    input:
        seq = expand("intermediate/ldhelm/ldhelm_in/{{spec}}/cpg_filt/{{spec}}.{scaff}.{{phase}}.seq", scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        pos = expand("intermediate/ldhelm/ldhelm_in/{{spec}}/cpg_filt/{{spec}}.{scaff}.{{phase}}.zero_base.pos", scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        anc_pos = expand('intermediate/ldhelm/anc_priors/{{spec}}/cpg_filt/{{spec}}.{scaff}.anc_priors.zero_pos.txt', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'])
    output:
        confs = "intermediate/ldhelm/conf/{spec}/{spec}.{phase}.conf"
    shell:
        "ldhelmet find_confs -w 50 -o {output.confs} {input.seq}"

## generate ldhelmet liklihood lookup table for all scaffs combined
rule ldhelm_lk:
    input:
        confs = "intermediate/ldhelm/conf/{spec}/cpg_filt/{spec}.{phase}.conf",
        seq = expand("intermediate/ldhelm/ldhelm_in/{{spec}}/cpg_filt/{{spec}}.{scaff}.{{phase}}.seq", scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        pos = expand("intermediate/ldhelm/ldhelm_in/{{spec}}/cpg_filt/{{spec}}.{scaff}.{{phase}}.zero_base.pos", scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        anc_pos = expand('intermediate/ldhelm/anc_priors/cpg_filt/{{spec}}/{{spec}}.{scaff}.anc_priors.zero_pos.txt', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'])
    output:
        lk = "intermediate/ldhelm/lk/{spec}/{spec}.{phase}.lk"
    shell:
        "ldhelmet table_gen --num_threads 20 -c {input.confs} -t 0.003 -r 0.0 0.1 10.0 1.0 100.0 -o {output.lk}"

## run ldhelmet pade for all scaffs combined
rule ldhelm_pade:
    input:
        confs = "intermediate/ldhelm/conf/{spec}/cpg_filt/{spec}.{phase}.conf",
        seq = expand("intermediate/ldhelm/ldhelm_in/{{spec}}/cpg_filt/{{spec}}.{scaff}.{{phase}}.seq", scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        pos = expand("intermediate/ldhelm/ldhelm_in/{{spec}}/cpg_filt/{{spec}}.{scaff}.{{phase}}.zero_base.pos", scaff = scaffolds.scaff_name[scaffolds['category']=='non_z']),
        anc_pos = expand('intermediate/ldhelm/anc_priors/{{spec}}/cpg_filt/{{spec}}.{scaff}.anc_priors.zero_pos.txt', scaff = scaffolds.scaff_name[scaffolds['category']=='non_z'])
    output:
        pade = "intermediate/ldhelm/pade/{spec}/{spec}.{phase}.pade"
    shell:
        "ldhelmet pade --num_threads 20 -c {input.confs} -t 0.003 -x 11 -o {output.pade}"

## run ldehlmet rjmcmc for each scaff separately for 5 iterations
rule ldhelm_rjmcmc:
    input:
        confs = "intermediate/ldhelm/conf/{spec}/cpg_filt/{spec}.{phase}.conf",
        seq = "intermediate/ldhelm/ldhelm_in/{spec}/cpg_filt/{spec}.{scaff}.{phase}.seq",
        pos = "intermediate/ldhelm/ldhelm_in/{spec}/cpg_filt/{spec}.{scaff}.{phase}.zero_base.pos",
        anc_pos = "intermediate/ldhelm/anc_priors/{spec}/no_cpg_filt/{spec}.{scaff}.anc_priors.zero_pos.no_cpg_filt.txt",
        mut_mat = "data/mut_mat.txt",
        lk = "intermediate/ldhelm/lk/{spec}/cpg_filt/{spec}.{phase}.lk",
        pade = "intermediate/ldhelm/pade/{spec}/cpg_filt/{spec}.{phase}.pade"
    output:
        rjmcmc = "intermediate/ldhelm/rjmcmc/{spec}/no_cpg_filt/{spec}.{scaff}.{phase}.{iter}.rjmcmc.post"
    shell:
        "ldhelmet rjmcmc --num_threads 10 -w 50 -l {input.lk} -p {input.pade} -b 10 --pos_file {input.pos} --snps_file {input.seq} -m {input.mut_mat} --burn_in 200000 -n 2000000 -a {input.anc_pos} -o {output.rjmcmc} --seed $RANDOM"

## Convert ldhelmet output to text 
rule post_to_text:
    input:
        rjmcmc = "intermediate/ldhelm/rjmcmc/{spec}/no_cpg_filt/{spec}.{scaff}.{phase}.{iter}.rjmcmc.post"
    output:
        "intermediate/ldhelm/rho/{spec}/no_cpg_filt/{spec}.{scaff}.{phase}.{iter}.rjmcmc.txt"
    shell:
        "ldhelmet post_to_text -m -p 0.025 -p 0.5 -p 0.975 -o {output} {input.rjmcmc}"


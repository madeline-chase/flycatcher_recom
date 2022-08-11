configfile: 'config.yaml'
localrules: all

rule all:
    input:
        vcfs = expand('data/vcf/all_sites.{scaff}.vcf.gz', scaff = config['scaffs'])

## Apply BQSR scores to bam files
rule recal_bam:
    input:
        bam = 'data/aln/{sample}.merged.marked.bam',
        table = 'data/aln/{sample}.recal_1.table',
        reference = config['reference']
    output:
        bam = 'data/aln/{sample}.merged.recal.bam'
    shell:
        """
        module load bioinfo-tools GATK/4.1.4.1
        gatk ApplyBQSR -R {input.reference} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam}
        """

## Run haplotype caller on bam files        
rule HC_from_bam:
    input:
        bam = 'data/aln/{sample}.merged.recal.bam',
        reference = config['reference']
    output:
        gvcf = 'data/gvcf/{sample}.{scaff}.g.vcf.gz'
    shell:
        """
        module load bioinfo-tools GATK/4.1.4.1
        gatk --java-options "-Xmx60g" HaplotypeCaller -R {input.reference} -I {input.bam} -O {output.gvcf} -ERC GVCF -L {wildcards.scaff}
        """

## Run haplotype caller on cram files
rule HC_from_cram:
    input:
        cram = 'data/aln/{sample}.merged.recal.cram',
        reference = config['reference']
    output:
        gvcf = 'data/gvcf/{sample}.{scaff}.g.vcf.gz'
    shell:
        """
        module load bioinfo-tools GATK/4.1.4.1
        gatk --java-options "-Xmx60g" HaplotypeCaller -R {input.reference} -I {input.cram} -O {output.gvcf} -ERC GVCF -L {wildcards.scaff}
        """

## Combine gvcfs
rule DB_import:
    input:
        map_file = 'data/mapfiles/{scaff}.mapfile.txt',
        gvcfs = expand('data/gvcf/{sample}.{{scaff}}.g.vcf.gz', sample = config['sample'])
    output:
        db_path = directory('data/gen_dbs/gen_db_{scaff}')
    shell:
        """
        module load bioinfo-tools GATK/4.1.4.1
        gatk --java-options "-Xmx40g -Xms40g" GenomicsDBImport --genomicsdb-workspace-path {output.db_path} -L {wildcards.scaff} --tmp-dir $SNIC_TMP --sample-name-map {input.map_file}
        """

## Run genotype gvcfs for all samples
rule gt_gvcf:
    input:
        db_path = 'data/gen_dbs/gen_db_{scaff}',
        reference = config['reference']
    output:
        vcf = 'data/vcf/all_sites.{scaff}.vcf.gz'
    shell:
        """
        module load bioinfo-tools GATK/4.1.4.1
        gatk --java-options "-Xmx30g" GenotypeGVCFs -R {input.reference} -V gendb://{input.db_path} -all-sites -L {wildcards.scaff} -O {output.vcf}
        """

"""
This rule takes any vcf file, compresses-it, then indexes it
"""
rule compress_vcf:
    input:
        "snpsift/GeneSets/{sample}.vcf"
    output:
        "snpsift/GeneSets/{sample}.vcf.gz"
    message:
        "Compressing and indexing {wildcards.sample}"
    threads:
        config["threads"]
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 90)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(8 * attempt, 15)
        )
    conda:
        "../envs/pbgzip.yaml"
    log:
        "logs/compress/{sample}.vcf.log"
    shell:
        "pbgzip -c {input} -n {threads} > {output} 2> {log}"


"""
This rule indexes the dbNSFP file if it does not exist
"""
rule index_dbnsfp:
    input:
        refs_dict["dbNSFP"]
    output:
        refs_dict["dbNSFP_tbi"]
    message:
        "Indexing dbNSFP file with tabix"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    params:
        extra = ""
    log:
        "logs/tabix/dbNSFP.log"
    wrapper:
        f"{swv}/bio/tabix"


"""
This rule indexes a compressed VCF file
"""
rule vcf_tabix:
    input:
        "snpsift/GeneSets/{sample}.vcf.gz"
    output:
        "snpsift/GeneSets/{sample}.vcf.gz.tbi"
    message:
        "Indexing annotated vcf for {wildcards.sample}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    params:
        extra = "-p vcf"
    log:
        "logs/tabix/{sample}.log"
    wrapper:
        f"{swv}/bio/tabix"

"""
This rule adds dbNSFP annotation on (annotated) variant calls
"""
rule snpsift_dbNSFP:
    input:
        calls = "snpeff/annotate/{sample}.vcf",
        dbNSFP = refs_dict["dbNSFP"],
        tbi = refs_dict["dbNSFP_tbi"]
    output:
        calls = temp("snpsift/dbNSFP/{sample}.vcf")
    message:
        "Adding dbNSFP annotation on {wildcards.sample} calls"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 90)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(8 * attempt, 15)
        )
    conda:
        "../envs/SnpSift.yaml"
    params:
        extra = config["params"].get("snpsift_dbNSFP_extra", "-v")
    log:
        "snpsift/logs/dbNSFP.{sample}.log"
    script:
        "../scripts/snpsift_dbNSFP.py"

"""
This rule adds gwas catalog annotation on (annotated) variant calls
"""
rule snpsift_GWASCat:
    input:
        calls = "snpsift/dbNSFP/{sample}.vcf",
        GWASCat = refs_dict["GWASCat"]
    output:
        calls = temp("snpsift/GWASCat/{sample}.vcf")
    message:
        "Adding GWAS Catalog annotations on {wildcards.sample} calls"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 90)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(8 * attempt, 15)
        )
    conda:
        "../envs/SnpSift.yaml"
    params:
        extra = config["params"].get("snpsift_GWASCat_extra", "-v")
    log:
        "snpsift/logs/GWASCat.{sample}.log"
    script:
        "../scripts/snpsift_GWASCat.py"

"""
This rule adds gene sets informations on (annotated) calls based on MSigDB
"""
rule snpsift_GeneSets:
    input:
        calls = "snpsift/GWASCat/{sample}.vcf",
        geneSets = refs_dict["GeneSets"]
    output:
        calls = temp("snpsift/GeneSets/{sample}.vcf")
    message:
        "Adding Genes Sets information on {wildcards.sample} based on MSigDB"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 90)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(8 * attempt, 15)
        )
    conda:
        "../envs/SnpSift.yaml"
    params:
        extra = config["params"].get("snpsift_GeneSets_extra", "-v")
    log:
        "snpsift/logs/GenesSets.{sample}.log"
    script:
        "../scripts/snpsift_GeneSets.py"

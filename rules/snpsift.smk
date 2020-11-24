"""
This rule add a variant type on a given VCF file based on the format
information and the ref/alt columns.
"""
rule snpsift_vartype:
    input:
        vcf = "snpeff/annotate/{sample}.vcf"
    output:
        vcf = "snpsift/vartype/{sample}.vcf"
    message:
        "Adding variant information on {wildcards.sample}"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 120)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(1024 * attempt, 8192)
        )
    params:
        extra = config["params"].get("snpsift_varType_extra", "-v")
    log:
        "logs/snpsift_vartype/{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/varType"


"""
This rule adds dbNSFP annotation on (annotated) variant calls
"""
rule snpsift_dbNSFP:
    input:
        vcf = "snpsift/vartype/{sample}.vcf",
        dbNSFP = refs_dict["dbNSFP"],
        tbi = refs_dict["dbNSFP_tbi"]
    output:
        vcf = temp("snpsift/dbNSFP/{sample}.vcf")
    message:
        "Adding dbNSFP annotation on {wildcards.sample} calls"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 120)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(1024 * attempt, 8192)
        )
    params:
        extra = config["params"].get("snpsift_dbNSFP_extra", "-v")
    log:
        "snpsift/logs/dbNSFP.{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/dbnsfp"


"""
This rule adds gwas catalog annotation on (annotated) variant calls
"""
rule snpsift_GWASCat:
    input:
        vcf = "snpsift/dbNSFP/{sample}.vcf",
        gwascat = refs_dict["GWASCat"]
    output:
        vcf = temp("snpsift/GWASCat/{sample}.vcf")
    message:
        "Adding GWAS Catalog annotations on {wildcards.sample} calls"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 120)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(1024 * attempt, 8192)
        )
    params:
        extra = config["params"].get("snpsift_GWASCat_extra", "-v")
    log:
        "snpsift/logs/GWASCat.{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/gwascat"


"""
This rule adds gene sets informations on (annotated) calls based on MSigDB
"""
rule snpsift_GeneSets:
    input:
        vcf = "snpsift/GWASCat/{sample}.vcf",
        gmt = refs_dict["GeneSets"]
    output:
        vcf = temp("snpsift/GeneSets/{sample}.vcf")
    message:
        "Adding Genes Sets information on {wildcards.sample} based on MSigDB"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 120)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(1024 * attempt, 8192)
        )
    params:
        extra = config["params"].get("snpsift_GeneSets_extra", "-v")
    log:
        "snpsift/logs/GenesSets.{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/genesets"

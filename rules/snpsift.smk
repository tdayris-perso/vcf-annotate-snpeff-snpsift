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
        time = (
            lambda wildcards, attempt: min(attempt * 30, 115)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 8192)
        )
    conda:
        "../envs/SnpSift.yaml"
    params:
        extra = config["params"].get("snpsift_dbNSFP_extra", "-v")
    log:
        "logs/snpsift/dbNSFP/{sample}.log"
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
        time = (
            lambda wildcards, attempt: min(attempt * 30, 115)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 8192)
        )
    conda:
        "../envs/SnpSift.yaml"
    params:
        extra = config["params"].get("snpsift_GWASCat_extra", "-v")
    log:
        "logs/snpsift/GWASCat/{sample}.log"
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
        time = (
            lambda wildcards, attempt: min(attempt * 30, 115)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 8192)
        )
    conda:
        "../envs/SnpSift.yaml"
    params:
        extra = config["params"].get("snpsift_GeneSets_extra", "-v")
    log:
        "logs/snpsift/GenesSets/{sample}.log"
    script:
        "../scripts/snpsift_GeneSets.py"


"""
This rule adds "SNP/MNP/INS/DEL/MIXED" in the INFO field.
It also adds "HOM/HET", but this last one works if there is only one sample
(otherwise it doesn't make any sense).
"""
rule snpsift_varType:
    input:
        calls = "snpsift/GeneSets/{sample}.vcf"
    output:
        calls = temp("snpsift/varType/{sample}.vcf")
    message:
        "Adding variant type in INFO field on {wildcards.sample}"
    threads:
        1
    resources:
        time = (
            lambda wildcards, attempt: min(attempt * 30, 115)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 8192)
        )
    conda:
        "../envs/SnpSift.yaml"
    log:
        "logs/snpsift/varType/{sample}.log"
    script:
        "../scripts/snpsift_varType.py"


"""
This rule annotates using PhastCons conservation scores.
"""
rule phastCons:
    input:
        call = "snpsift/varType/{sample}.vcf",
        database = config["ref"]["phastCons"]
    output:
        call = temp("snpsift/phastCons/{sample}.vcf")
    message:
        "Annotating {wildcards.sample} with phastCons database"
    threads:
        1
    resources:
        time = (
            lambda wildcards, attempt: min(attempt * 30, 115)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 8192)
        )
    conda:
        "../envs/SnpSift.yaml"
    log:
        "logs/snpsift/varType/{sample}.log"
    script:
        "../scripts/snpsift_varType.py"


"""
Annotate using fields from another VCF file
(e.g. dbSnp, 1000 Genomes projects, ClinVar, ExAC, etc.)
"""
rule snpsift_dbsnp:
    input:
        call = "snpeff/annotate/{sample}.vcf",
        database = config["ref"]["dbSNP"],
        database_index = f"config['ref']['dbSNP'].tbi"
    output:
        call = "snpsift/dbsnp/{sample}.vcf"
    message:
        "Annotating {sample} with dbsnp"
    threads:
        1
    resources:
        time = (
            lambda wildcards, attempt: min(attempt * 30, 115)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 8192)
        )
    conda:
        "../envs/SnpSift.yaml"
    params:
        extra = config["params"].get(
            "snpsift_dbSNP_extra",
            "-a -tabix -noDownload"
        )
    log:
        "logs/snpsift/dbsnp/{sample}.log"
    script:
        "../scripts/snpsift_dbsnp.py"

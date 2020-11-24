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
        call = "snpsift/vartype/{sample}.vcf",
        dbNSFP = refs_dict["dbNSFP"],
        tbi = refs_dict["dbNSFP_tbi"]
    output:
        call = temp("snpsift/dbNSFP/{sample}.vcf")
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
        "logs/snpsift/dbNSFP/{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/dbnsfp"


"""
This rule adds gwas catalog annotation on (annotated) variant calls
"""
rule snpsift_GWASCat:
    input:
        call = "snpsift/dbNSFP/{sample}.vcf",
        gwascat = refs_dict["GWASCat"]
    output:
        call = temp("snpsift/GWASCat/{sample}.vcf")
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
        "logs/snpsift/GWASCat/{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/gwascat"


"""
This rule adds gene sets informations on (annotated) calls based on MSigDB
"""
rule snpsift_GeneSets:
    input:
        call = "snpsift/GWASCat/{sample}.vcf",
        gmt = refs_dict["GeneSets"]
    output:
        call = temp("snpsift/GeneSets/{sample}.vcf")
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
        "logs/snpsift/GenesSets/{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/genesets"


"""
This rule adds 1000g annotation based on its VCF.
"""
rule snpsift_1000genomes:
    input:
        call = "snpsift/GeneSets/{sample}.vcf",
        vcf = refs_dict["1000genomes"]
    output:
        call = temp("snpsift/1000genomes/{sample}.vcf")
    message:
        "Annotating with 1000 genomes variants"
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
        extra = config["params"].get("snpsift_annotate_extra", "-v")
    log:
        "logs/snpsift/1000g/{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/annotate"

"""
This rule adds clinvar annotation based on its VCF.
"""
rule snpsift_clinvar:
    input:
        call = "snpsift/1000genomes/{sample}.vcf",
        vcf = refs_dict["clinvar"]
    output:
        call = temp("snpsift/clinvar/{sample}.vcf")
    message:
        "Annotating with ClinVar database"
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
        extra = config["params"].get("snpsift_annotate_extra", "-v")
    log:
        "logs/snpsift/clinvar/{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/annotate"

"""
This rule adds cosmic annotation based on its VCF.
"""
rule snpsift_cosmic:
    input:
        call = "snpsift/clinvar/{sample}.vcf",
        vcf = refs_dict["cosmic"]
    output:
        call = temp("snpsift/Cosmic/{sample}.vcf")
    message:
        "Annotating with Cosmic database"
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
        extra = config["params"].get("snpsift_annotate_extra", "-v")
    log:
        "logs/snpsift/cosmic/{sample}.log"
    wrapper:
        f"{git}/bio/snpsift/annotate"

"""
This rule calls snpeff to annotate from a VCF file.
More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/snpeff.html
"""
rule snpeff:
    input:
        vcf = (
            lambda wildcards: get_vcf_w(wildcards.sample)
        )
    output:
        calls = "snpeff/annotate/{sample}.vcf",
        stats = report(
            "snpeff/report/{sample}.html",
            caption="../report/snpeff_report.rst",
            category="Quality"
        ),
        csvstats = "snpeff/stats/{sample}.csv"
    message:
        "Annotating {wildcards.sample} with SnpEff"
    threads:
        1
    resources:
        time = (
            lambda wildcards, attempt: min(attempt * 30, 90)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8, 15)
        )
    params:
        extra = config["params"].get("snpeff_extra", "-v"),
        reference = config["params"].get("organism", "GRCh38.86")
    log:
        "snpeff/logs/{sample}.log"
    wrapper:
        f"{swv}/bio/snpeff"

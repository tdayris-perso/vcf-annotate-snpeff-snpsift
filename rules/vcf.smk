"""
This rule takes any vcf file, compresses-it, then indexes it
"""
rule compress_vcf:
    input:
        "snpsift/GeneSets/{sample}.vcf"
    output:
        report(
            "snpsift/GeneSets/{sample}.vcf.gz",
            caption="../report/vcf_annotated.rst",
            category="Calls"
        )
    message:
        "Compressing and indexing {wildcards.sample}"
    threads:
        1
    resources:
        time_min = (
            lambda wildcars, attempt: min(30 * attempt, 90)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(8 * attempt, 15)
        )
    wrapper:
        f"{swv}/bio/vcf/compress"

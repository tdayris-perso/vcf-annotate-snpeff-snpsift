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
        time = (
            lambda wildcards, attempt: min(attempt * 30, 115)
        ),
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 128, 512)
        )
    log:
        "logs/compress/{sample}.vcf.log"
    shell:
        "gzip -c {input} > {output} 2> {log}"

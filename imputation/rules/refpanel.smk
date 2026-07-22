#creates the reference panel for imputation
__author__ = "Evelyn Todd"
__copyright__ = "Copyright (c) 2026 Evelyn Todd"
__license__ = "MIT"

# ---- final outputs ----
rule all_refpanel:
    input:
        expand("refpanel/ref_panel_sites_chr{chr}.vcf.gz", chr=config["chromosomes"]),
        expand("refpanel/ref_panel_sites_chr{chr}.tsv.gz", chr=config["chromosomes"]),
        expand("refpanel/chunks/chunks.chr{chr}.txt", chr=config["chromosomes"])      

# ---- extract sites from ref panel ----
rule extract_sites:
    input:
        canidvcf = lambda wc: config["canidvcf"].format(chr=wc.chr)
    output:
        sitesvcf = "refpanel/ref_panel_sites_chr{chr}.vcf.gz",
        sitesvcfidx = "refpanel/ref_panel_sites_chr{chr}.vcf.gz.csi",
        sitestsv = "refpanel/ref_panel_sites_chr{chr}.tsv.gz"
    resources:
        mem_mb=50000,
        runtime=4*60
    threads: 4
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps {input.canidvcf} --threads {threads} -Oz -o {output.sitesvcf}
        bcftools index -f {output.sitesvcf}
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sitesvcf} | bgzip -c > {output.sitestsv}
        tabix -s1 -b2 -e2 {output.sitestsv}
        """

# ---- create chunks for ref panel ----
rule glimpse_chunk:
    input:
        sitesvcf = "refpanel/ref_panel_sites_chr{chr}.vcf.gz",
    output:
        chunks = "refpanel/chunks/chunks.chr{chr}.txt",
    resources:
        mem_mb=20000,
        runtime=1*60
    conda:
        "envs/glimpse.yaml"
    shell:
        """
        GLIMPSE_chunk --input {input.sitesvcf} --region chr{wildcards.chr} --window-size 2000000 --buffer-size 200000 --output {output.chunks}
        """
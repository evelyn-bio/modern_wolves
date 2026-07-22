# downsample bam files and test imputation
__author__ = "Evelyn Todd"
__copyright__ = "Copyright (c) 2026 Evelyn Todd"
__license__ = "MIT"


# Read file with filepath, sample name and coverage
bam_paths = {}
coverage_vals = {}
samples = []
with open(config["downsamplebams"]) as f:
    for line in f:
        line = line.strip()
        if not line:  # skip empty lines
            continue
        bam, sample, cov = line.split("\t")
        bam_paths[sample] = bam
        coverage_vals[sample] = cov
        samples.append(sample)
fractions = [str(f) for f in config["downsamplefractions"]]

# ---- final outputs ----
rule all_downsample:
    input:
        expand(["downsample/bams/{sample}_{fraction}.bam","downsample/bams/{sample}_{fraction}.coverage"],sample=samples,fraction=fractions),
        expand(["downsample/variants/{sample}_{fraction}_chr1.vcf.gz"],sample=samples,fraction=fractions),
        expand("downsample/ligate/{sample}_{fraction}_chr1_ligate.bcf",sample=samples,fraction=fractions),
        expand(["downsample/merge/{sample}_GP{gp}.vcf.gz"],sample=samples,gp=config["gp_thresholds"]),
        expand(["downsample/stats/{sample}_GP{gp}.imiss","downsample/stats/{sample}_GP{gp}.genome",
        "downsample/stats/{sample}_GP{gp}.het","downsample/stats/{sample}_GP{gp}.hom.indiv"],sample=samples,gp=config["gp_thresholds"])


# ---- downsample bams ----
rule downsample_bam:
    input:
        bam=lambda wc: bam_paths[wc.sample]
    output:
        bam="downsample/bams/{sample}_{fraction}.bam",
        cov="downsample/bams/{sample}_{fraction}.coverage"
    params:
        prop=lambda wc: float(wc.fraction) / float(coverage_vals[wc.sample])
    resources:
      mem_mb=50000,
      runtime=2*60
    shell:
        """
        [ -f {output.bam} ] && rm {output.bam}
        picard DownsampleSam -I {input.bam} -O {output.bam} -P {params.prop} -CREATE_INDEX True -STRATEGY Chained
        paleomix coverage --overwrite-output {output.bam} {output.cov}
        tail -n +22 {output.cov} | head -1 | awk '{{print "{wildcards.sample}", $NF}}' >> downsample/bams/coverage.txt
        """

# ---- variant calling ----
rule call_variants:
    input:
        sitesvcf = "refpanel/ref_panel_sites_chr1.vcf.gz",
        sitestsv = "refpanel/ref_panel_sites_chr1.tsv.gz",
        bamlist = "downsample/bams/{sample}_{fraction}.bam",
        ref = config["ref_genome"]
    output:
        targetvcf = "downsample/variants/{sample}_{fraction}_chr1.vcf.gz",
        targetidx = "downsample/variants/{sample}_{fraction}_chr1.vcf.gz.csi"
    resources:
        mem_mb=40000,
        runtime=3*60
    threads: 4
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -I -E -a 'FORMAT/DP' -T {input.sitesvcf} -r chr1 {input.bamlist} -Ou | \
        bcftools call -Aim -C alleles -T {input.sitestsv} -Oz -o {output.targetvcf}
        bcftools index -f {output.targetvcf}
        """

# ---- phase each separately----
rule glimpse_phase:     
    input:
        targetvcf = "downsample/variants/{sample}_{fraction}_chr1.vcf.gz",
        canidvcf="/projects/psg/people/pkb156/AW/canidref/ref-panel_chr1_sample-snp_filltags_filter.phased.vcf.gz",
        map = "/projects/psg/people/pkb156/AW_old/downsample/glimpse/maps/chr1_average_canFam3.1.txt",
        chunks="refpanel/chunks/chunks.chr1.txt"
    output:
        imputed = "downsample/phase/{sample}_{fraction}_chr1.imputed.00.bcf"
    params:    
        prefix = "downsample/phase/{sample}_{fraction}_chr1.imputed"
    resources:
      mem_mb=50000,
      runtime=5*60
    threads: 2
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT={params.prefix}.${{ID}}.bcf
            GLIMPSE_phase \
            --input {input.targetvcf} \
            --reference {input.canidvcf} \
            --map {input.map} \
            --input-region ${{IRG}} \
            --output-region ${{ORG}} --output ${{OUT}} \
            --thread {threads}
            bcftools index -f ${{OUT}}
        done < {input.chunks}
        touch {output.imputed}
        """
        
# ----list of chunks to ligate----    
rule ligate_list:
    input:
        chunks = "refpanel/chunks/chunks.chr1.txt",
        imputed = "downsample/phase/{sample}_{fraction}_chr1.imputed.00.bcf"
    output:
        ligated_list = "downsample/ligate/{sample}_{fraction}_chr1_ligated_list.txt"
    resources:
      mem_mb=5000,
      runtime=10
    params:
        prefix = "downsample/phase/{sample}_{fraction}_chr1.imputed"
    shell:
        """
        > {output.ligated_list}
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            ls {params.prefix}.${{ID}}.bcf >> {output.ligated_list}
        done < {input.chunks}
        """

# ----ligate into chr----    
rule ligate_chr:
    input:
        ligated_list = "downsample/ligate/{sample}_{fraction}_chr1_ligated_list.txt"
    output:
        ligated_bcf = "downsample/ligate/{sample}_{fraction}_chr1_ligate.bcf"
    resources:
      mem_mb=10000,
      runtime=60
    threads: 8
    shell:
        """
        GLIMPSE_ligate \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads}
        bcftools index -f {output.ligated_bcf}
        """

# ----filter by GP----    
rule GP_filter:
    input:
        ligated_bcf = "downsample/ligate/{sample}_{fraction}_chr1_ligate.bcf"
    output:
        vcf = "downsample/GP/{sample}_{fraction}_GP{gp}.vcf.gz"
    resources:
        mem_mb=50000,
        runtime=2*60
    threads: 1
    shell:
        """
        # Filter by GP
        bcftools filter -Oz {input.ligated_bcf} -i 'FORMAT/GP >= {wildcards.gp}' -o {output.vcf}
        bcftools index -f {output.vcf}
        """
        
# ---- merge filtered fractions into one VCF per sample & GP ----
rule merge_files:
    input:
        gp_vcf = lambda wc: expand(
            "downsample/GP/{sample}_{fraction}_GP{gp}.vcf.gz",
            sample=[wc.sample],
            fraction=fractions,
            gp=[wc.gp]
        )
    output:
        merged_vcf = "downsample/merge/{sample}_GP{gp}.vcf.gz"
    threads: 1
    resources:
        mem_mb=40000,
        runtime=2*60
    shell:
        """
        bcftools merge --force-samples -Oz {input.gp_vcf} -o {output.merged_vcf}
        bcftools index -f {output.merged_vcf}
        """

# ---- missing and similarity checks ----
rule plink_stats:
    input:
        merged_vcf = "downsample/merge/{sample}_GP{gp}.vcf.gz"
    output:
        miss = "downsample/stats/{sample}_GP{gp}.imiss",
        related = "downsample/stats/{sample}_GP{gp}.genome",
        hom = "downsample/stats/{sample}_GP{gp}.hom.indiv"
    params:
        filename = "downsample/stats/{sample}_GP{gp}"
    threads: 1
    resources:
        mem_mb=40000,
        runtime=2*60
    shell:
        """
        plink --vcf {input.merged_vcf} --dog --allow-extra-chr --missing --out {params.filename}
        plink --vcf {input.merged_vcf} --dog --allow-extra-chr --genome --out {params.filename}
        plink --vcf {input.merged_vcf} --dog --allow-extra-chr --homozyg --out {params.filename}
        """
        
       
# ---- heterozygosity per sample ----        
rule het_per_sample:
    input:
        vcf = "downsample/merge/{sample}_GP{gp}.vcf.gz"
    output:
        het="downsample/stats/{sample}_GP{gp}.het"
    threads: 1
    resources:
        mem_mb = 40000,
        runtime = 120
    shell:
        """
        bcftools stats -s - {input.vcf} > {output.het}
        """


# joint impute all samples in target panel
__author__ = "Evelyn Todd"
__copyright__ = "Copyright (c) 2026 Evelyn Todd"
__license__ = "MIT"

# ---- final outputs ----
rule all_imputepanel:
    input:
        expand(["impute/variants/targetpanel_chr{chr}.vcf.gz", "impute/variants/targetpanel_chr{chr}.vcf.gz.csi"],chr=config["chromosomes"]),
        expand(["impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids.bcf"],chr=config["chromosomes"]),
        expand(["impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore08.bcf"],chr=config["chromosomes"]),
        "impute/subset/targetpanel_ligate_reheader_wolvescanids_infoscore08_counts.txt",
        expand(["impute/subset/GP0.90_90_chr{chr}.txt"],chr=config["chromosomes"]),
        expand(["impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter.vcf.gz"],chr=config["chromosomes"]),
        expand(["impute/refpanel/refpanel_samples_chr{chr}_GPfilter.vcf.gz"],chr=config["chromosomes"]),
        expand(["impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter_phase.vcf.gz"],chr=config["chromosomes"]),
        expand(["impute/combine/combined_chr{chr}.vcf.gz"],chr=config["chromosomes"]),
        expand(["impute/combine/combined_chr{chr}_targetpanel.vcf.gz"],chr=config["chromosomes"]),
        "impute/combine/combined_allchr.vcf.gz",
        expand("impute/variants/targetpanel_chr{chr}.vcf.gz",chr=config["chromosomes"]),
        expand(["impute/combine/combined_chr{chr}_infoscore.vcf.gz"],chr=config["chromosomes"]),
        "impute/combine/combined_allchr_infoscore.vcf.gz"         
        

#split -l 200 data/bamfiles.list data/bam_chunk_ #split bam files into chunks of 200 to run separately
import glob
CHUNKS = [os.path.basename(x).replace("data/", "") 
          for x in glob.glob("data/bam_chunk_*")]

# ---- variant calling ----       
rule variants_impute_chunk:
    input:
        sitesvcf = "refpanel/ref_panel_sites_chr{chr}.vcf.gz",
        sitestsv = "refpanel/ref_panel_sites_chr{chr}.tsv.gz",
        bamlist = "data/{chunk}",
        ref = config["ref_genome"]
    output:
        vcf=temp("impute/tmp/chr{chr}.{chunk}.bcf"),
        vcfidx=temp("impute/tmp/chr{chr}.{chunk}.bcf.csi")
    threads: 4
    resources:
        mem_mb=150000,
        runtime=25*60
    shell:
        """
        bcftools mpileup --ignore-RG --threads {threads} \
        -f {input.ref} -I -E -q 30 -Q 20 \
        -a FORMAT/DP \
        -T {input.sitesvcf} \
        -r chr{wildcards.chr} \
        -b {input.bamlist} -Ou | \
        bcftools call -Aim -C alleles \
        -T {input.sitestsv} -Ob -o {output.vcf}
        bcftools index -f {output.vcf}
        """
rule merge_chunks:
    input:
        vcf=lambda wc: expand("impute/tmp/chr{chr}.{chunk}.bcf",chr=wc.chr,chunk=CHUNKS),
        vcfidx=lambda wc: expand("impute/tmp/chr{chr}.{chunk}.bcf.csi",chr=wc.chr,chunk=CHUNKS)
    output:
        vcf="impute/variants/targetpanel_chr{chr}.vcf.gz",
        vcfidx="impute/variants/targetpanel_chr{chr}.vcf.gz.csi"
    threads: 4
    resources:
        mem_mb=80000,
        runtime=8*60
    shell:
        """
        bcftools merge --force-samples --threads {threads} {input.vcf} -Oz -o {output.vcf}
        bcftools index -f {output.vcf}
        """

# ---- phase by chr ----
rule impute_phase:     
    input:
        targetvcf = "impute/variants/targetpanel_chr{chr}.vcf.gz",
        canidvcf = lambda wc: config["canidvcf"].format(chr=wc.chr),
        map =  lambda wc: config["recmap"].format(chr=wc.chr),
        chunks="refpanel/chunks/chunks.chr{chr}.txt"
    output:
        imputed = "impute/phase/targetpanel_chr{chr}.imputed.00.bcf"
    params:    
        prefix = "impute/phase/targetpanel_chr{chr}.imputed"
    resources:
      mem_mb=100000,
      runtime=20*60
    threads: 2
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        OUT={params.prefix}.${{ID}}.bcf
    
        if [ ! -f "${{OUT}}" ]; then
          echo "Processing chunk ${{ID}}..."
          GLIMPSE_phase \
          --input {input.targetvcf} \
          --reference {input.canidvcf} \
          --map {input.map} \
          --input-region ${{IRG}} \
          --output-region ${{ORG}} --output ${{OUT}} \
          --thread {threads}
        else
          echo "Chunk ${{ID}} already exists, skipping..."
        fi
        done < {input.chunks}
        touch {output.imputed}
        """

# ---- index the files----
rule impute_index:     
    input:
        imputed = "impute/phase/targetpanel_chr{chr}.imputed.00.bcf"
    output:
        imputedidx = "impute/phase/targetpanel_chr{chr}.imputed.00.bcf.csi"
    resources:
      mem_mb=10000,
      runtime=2*60
    shell:
        """
        bcftools index -f {input.imputed}
        """

# ----list of chunks to ligate----    
rule impute_ligate_list:
    input:
        chunks = "refpanel/chunks/chunks.chr{chr}.txt",
        imputed = "impute/phase/targetpanel_chr{chr}.imputed.00.bcf"
    output:
        ligated_list = "impute/ligate/targetpanel_chr{chr}_ligated_list.txt"
    resources:
      mem_mb=5000,
      runtime=10
    params:
        prefix = "impute/phase/targetpanel_chr{chr}.imputed"
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
rule impute_ligate:
    input:
        ligated_list = "impute/ligate/targetpanel_chr{chr}_ligated_list.txt"
    output:
        ligated_bcf = "impute/phase/targetpanel_chr{chr}_ligate.bcf",
        ligated_idx = "impute/phase/targetpanel_chr{chr}_ligate.bcf.csi"
    resources:
      mem_mb=30000,
      runtime=5*60
    threads: 2
    shell:
        """
        GLIMPSE_ligate \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads}
        bcftools index -f {output.ligated_bcf}
        """
        
# ----reheader with sample names----    
rule reheader:
    input:
        bcf = "impute/phase/targetpanel_chr{chr}_ligate.bcf",
        sample_list = "data/bamfile_names.txt"
    output:
        bcf = "impute/phase/targetpanel_chr{chr}_ligate_reheader.bcf",
        idx = "impute/phase/targetpanel_chr{chr}_ligate_reheader.bcf.csi"
    resources:
      mem_mb=30000,
      runtime=1*60
    threads: 1
    shell:
        """
        bcftools reheader  --samples {input.sample_list} {input.bcf} -o {output.bcf}
        bcftools index -f {output.bcf}
        """

# ----subset samples of interest----    
rule subset_samples:
    input:
        bcf = "impute/phase/targetpanel_chr{chr}_ligate_reheader.bcf",
        sample_list = "data/wolvesandcanids.txt"
    output:
        bcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids.bcf",
        idx = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids.bcf.csi"
    resources:
      mem_mb=60000,
      runtime=4*60
    threads: 4
    shell:
        """
        bcftools view -Ob --threads {threads} --force-samples -S {input.sample_list} {input.bcf} -o {output.bcf}
        bcftools index -f {output.bcf}
        """

# ----filter for an infoscore of 0.8 ----    
rule infoscore_filter:
    input:
        bcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids.bcf"
    output:
        bcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore08.bcf",
        idx = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore08.bcf.csi",
        stats = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore08.stats"
    resources:
      mem_mb=60000,
      runtime=3*60
    threads: 4
    shell:
        """
        bcftools view -Ob --threads {threads} -i 'INFO>=0.8' {input.bcf} -o {output.bcf}
        bcftools index -f {output.bcf}
        bcftools stats {output.bcf} > {output.stats}
        """

# ----number of variants after infoscore filtering ----    
rule variant_counting:
    input:
        stats = expand("impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore08.stats",chr=config["chromosomes"])
    output:
        counts = "impute/subset/targetpanel_ligate_reheader_wolvescanids_infoscore08_counts.txt"
    resources:
        mem_mb=1000,
        runtime=5
    threads: 1
    shell:
        """
        for f in {input.stats}; do
            sed -n '24p' "$f"
        done > {output.counts}
        """

# ----filter rows by GP ----    
rule GP_filter_row:
    input:
        bcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids.bcf",
        pyscript = "scripts/GPfilter.py"
    output:
        bedfilter = "impute/subset/GP0.90_90_chr{chr}.txt"
    resources:
        mem_mb=8000,
        runtime=90
    threads: 1
    conda:
        "envs/gpfilter.yaml"
    shell:
        """
        module purge
        python {input.pyscript} {input.bcf} {output.bedfilter} 0.9 0.90
        """

# ----filter rows by GP ----    
rule remove_lowquality:
    input:
        bcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids.bcf",
        bedfilter = "impute/subset/GP0.90_90_chr{chr}.txt"
    output:
        filtervcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter.vcf.gz",
        filtervcfidx = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter.vcf.gz.csi",
    resources:
        mem_mb=50000,
        runtime=3*60
    threads: 1
    shell:
        """
        bcftools view -Oz -R {input.bedfilter} {input.bcf} -o {output.filtervcf}
        bcftools index -f {output.filtervcf}
        """


# ----remove samples of interest from reference panel ----    
rule subset_refpanel:
    input:
        refpanel = lambda wc: config["canidvcf"].format(chr=wc.chr),
        bedfilter = "impute/subset/GP0.90_90_chr{chr}.txt",
        samples = "data/refpanelsamples.txt",
        samplenames = "data/refpanelsamplesrename.txt"
    output:
        filtervcf = "impute/refpanel/refpanel_samples_chr{chr}_GPfilter.vcf.gz",
        filtervcfidx = "impute/refpanel/refpanel_samples_chr{chr}_GPfilter.vcf.gz.csi",
    resources:
        mem_mb=50000,
        runtime=3*60
    threads: 1
    shell:
        """
        bcftools view -Oz --force-samples -R {input.bedfilter} -S {input.samples} {input.refpanel} -o impute/refpanel/refpanel_samples_chr{wildcards.chr}_GPfilter_tmp.vcf.gz
        bcftools reheader  --samples {input.samplenames} impute/refpanel/refpanel_samples_chr{wildcards.chr}_GPfilter_tmp.vcf.gz -o {output.filtervcf}
        bcftools index -f {output.filtervcf}
        """

# ---- haplotype phasing of target panel ----  
rule phase_variants:
    input:
        filtervcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter.vcf.gz",
        bedfilter = "impute/subset/GP0.90_90_chr{chr}.txt"
    output:
        phasevcf = "impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter_phase.vcf.gz",
        phasevcfidx = "impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter_phase.vcf.gz.csi",
    resources:
        mem_mb=20000,
        runtime=60
    threads: 2
    conda:
        "envs/glimpse.yaml"
    shell:
        """
        GLIMPSE_sample --input {input.filtervcf} --solve --output {output.phasevcf} --thread {threads}
        bcftools index -f {output.phasevcf}
        """
 
 # ---- merge target and ref panels together ----  
rule combine_vcfs:
    input:
        targetvcf = "impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter_phase.vcf.gz",
        refvcf = "impute/refpanel/refpanel_samples_chr{chr}_GPfilter.vcf.gz"
    output:
        combinedvcf = "impute/combine/combined_chr{chr}.vcf.gz",
        combinedvcfidx = "impute/combine/combined_chr{chr}.vcf.gz.csi",
    resources:
        mem_mb=50000,
        runtime=3*60
    threads: 4
    shell:
        """
        bcftools merge -Oz --force-samples --threads {threads} {input.targetvcf} {input.refvcf} -o {output.combinedvcf}
        bcftools index -f {output.combinedvcf}
        """

# ---- concatenate all chromosomes together ----
rule concat_combine:
    input:
        vcfs = expand("impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_GPfilter_phase.vcf.gz", chr=config["chromosomes"])
    output:
        vcf = "impute/combine/combined_allchr_targetpanel.vcf.gz",
        idx = "impute/combine/combined_allchr_targetpanel.vcf.gz.csi",
        stats = "impute/combine/combined_allchr_targetpanel.vcf.gz.stats"
    resources:
        mem_mb=80000,
        runtime=3*60
    threads: 4
    shell:
        """
        bcftools concat --threads {threads} -Oz {input.vcfs} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.stats}
        """

# ----remove samples of interest from reference panel by infoscore ----    
rule subset_refpanel_infofilter:
    input:
        refpanel = lambda wc: config["canidvcf"].format(chr=wc.chr),
        targetpanel = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore08.bcf",
        samples = "data/refpanelsamples.txt",
        samplenames = "data/refpanelsamplesrename.txt"
    output:
        filtervcf = "impute/refpanel/refpanel_samples_chr{chr}_infoscore.vcf.gz",
        filtervcfidx = "impute/refpanel/refpanel_samples_chr{chr}_infoscore.vcf.gz.csi",
        positions = "impute/subset/targetpanel_chr{chr}.txt"
    resources:
        mem_mb=50000,
        runtime=3*60
    threads: 1
    shell:
        """
        bcftools query -f'%CHROM\t%POS\n' {input.targetpanel} | bgzip -c > {output.positions}
        bcftools view -Oz --force-samples -T {output.positions} -S {input.samples} {input.refpanel} -o impute/refpanel/refpanel_samples_chr{wildcards.chr}_GPfilter_tmp.vcf.gz
        bcftools reheader  --samples {input.samplenames} impute/refpanel/refpanel_samples_chr{wildcards.chr}_GPfilter_tmp.vcf.gz -o {output.filtervcf}
        bcftools index -f {output.filtervcf}
        """

# ---- haplotype phasing of target panel ----  
rule phase_variants_infoscore:
    input:
        filtervcf = "impute/subset/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore08.bcf"
    output:
        phasevcf = "impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore_phase.vcf.gz",
        phasevcfidx = "impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore_phase.vcf.gz.csi",
    resources:
        mem_mb=20000,
        runtime=60
    threads: 2
    conda:
        "envs/glimpse.yaml"
    shell:
        """
        GLIMPSE_sample --input {input.filtervcf} --solve --output {output.phasevcf} --thread {threads}
        bcftools index -f {output.phasevcf}
        """
 
 # ---- merge target and ref panels together ----  
rule combine_vcfs_infoscore:
    input:
        targetvcf = "impute/sample/targetpanel_chr{chr}_ligate_reheader_wolvescanids_infoscore_phase.vcf.gz",
        refvcf = "impute/refpanel/refpanel_samples_chr{chr}_infoscore.vcf.gz"
    output:
        combinedvcf = "impute/combine/combined_chr{chr}_infoscore.vcf.gz",
        combinedvcfidx = "impute/combine/combined_chr{chr}_infoscore.vcf.gz.csi",
    resources:
        mem_mb=50000,
        runtime=3*60
    threads: 4
    shell:
        """
        bcftools merge -Oz --force-samples --threads {threads} {input.targetvcf} {input.refvcf} -o {output.combinedvcf}
        bcftools index -f {output.combinedvcf}
        """

# ---- concatenate all chromosomes together ----
rule concat_combine_infoscore:
    input:
        vcfs = expand("impute/combine/combined_chr{chr}_infoscore.vcf.gz", chr=config["chromosomes"])
    output:
        vcf = "impute/combine/combined_allchr_infoscore.vcf.gz",
        idx = "impute/combine/combined_allchr_infoscore.vcf.gz.csi",
        stats = "impute/combine/combined_allchr_infoscore.vcf.gz.stats"
    resources:
        mem_mb=80000,
        runtime=3*60
    threads: 4
    shell:
        """
        bcftools concat --threads {threads} -Oz {input.vcfs} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.stats}
        """

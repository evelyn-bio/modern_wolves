#snakemake  --use-conda --use-env --workflow-profile slurm
workdir: "/projects/psg/people/pkb156/MW"
configfile: "config.yaml"

#load modules
shell.prefix("module load " + " ".join(config["modules"]) + " || true; ")


rule all:
    input:
        "popgen/phased.european.info08.name.norelated.maf01.vcf.gz",
        "popgen/phased.european.info08.name.nodups.maf01.vcf.gz",
        "popgen/phased.all.info08.name.norelated.maf01.vcf.gz",
        "popgen/phased.all.info08.name.nodups.maf01.vcf.gz",
        "smartpca/phased.european.info08.name.nodups.maf01.evec",
        "smartpca/phased.european.info08.name.nodups.maf01.eval",
        "popgen/phased.european.info08.name.norelated.maf01.wolves.vcf.gz",
        "admixture/phased.european.info08.name.norelated.maf01.prune.ped",
        "admixture/phased.european.info08.name.norelated.maf01.prune.map",
        expand("admixture/phased.european.info08.name.norelated.maf01.prune_{K}.out", K=config["K_range"]), 
        expand("haplonet/europeanall/haplonet.admixture.k{K}.f.npy", K=config["K_range"]),
        expand("haplonet/europeanall/haplonet.admixture.k{K}.q", K=config["K_range"]),
        "haplonet/europeanall/names.list",
        expand("PH/angsd_genotypes/{chr}_ancient_PH.bed", chr=config["chromosomes"]),
        expand("PH/angsd_genotypes/{chr}_ancient_PH.bim", chr=config["chromosomes"]),
        expand("PH/angsd_genotypes/{chr}_ancient_PH.fam", chr=config["chromosomes"]),
        "Dstats/dog/Dstat_dog_NWIB.rds",
        "Dstats/dog/Dstat_dog_SCAN.rds",
        "Dstats/wild/Dstat_wild_NWIB.rds",
        "Dstats/wild/Dstat_wild_SCAN.rds",
        "Dstats/ancient/Dstat_ancient_SCAN.rds",
        "Dstats/ancient/Dstat_ancient_NWIB.rds",
        "Dstats/dogdog/Dstat_dogdog.rds",
        "Dstats/asian/Dstat_asian_SCAN.rds",
        "Dstats/asian/Dstat_asian_NWIB.rds",
        "Dstats/asianregions/Dstat_asianregions.rds",
        "Dstats/ancientdog/Dstat_ancientdog.rds",
        "qpadm/f4ratio_all.RData",
        expand("orientagraph/og_rep/allpop_m{m}_rep{r}.treeout.gz", m=config["migrations"], r=config["replicates"]),
        expand("orientagraph/og_log/log_allpop_m{m}_rep{r}.log", m=config["migrations"], r=config["replicates"]),
        expand("orientagraph/og_rep2/allpop_m{m}_rep{r}.treeout.gz", m=config["migrations"], r=config["replicates"]),
        expand("orientagraph/og_log2/log_allpop_m{m}_rep{r}.log", m=config["migrations"], r=config["replicates"])
      
# ---- RELATED: identify related samples in king ----
rule king_related:
    input: 
        vcf=config["vcf"]
    output: 
        dups="king/removedups.txt",
        related="king/removerelated.txt",
        bed="popgen/phased.all.info08.name.bed",
        bim="popgen/phased.all.info08.name.bim",
        fam="popgen/phased.all.info08.name.fam"
    resources:
        mem_mb=80000,
        runtime=3*60
    shell:
        """
        plink --bcf {input.vcf} --double-id --dog --make-bed --out popgen/phased.all.info08.name
        king -b popgen/phased.all.info08.name.bed --related --sexchr 39 --prefix king/phased.all.info08.name --degree 2
        awk '$14=="Dup/MZ"{{print $4}}' king/phased.all.info08.name.kin0 > {output.dups}
        awk 'NR>1 && $10>0.177{{if(!seen[$3]++) print $3}}' king/phased.all.info08.name.kin0 > {output.related}
        """

# ---- FILE FILTER: remove ancient canids ----
rule remove_ancient:
    input:
        vcf=config["vcf"],
        wolves=config["wolves_list"]
    output: 
        vcf="popgen/phased.all.info08.name.noancient.vcf.gz",
        vcfidx="popgen/phased.all.info08.name.noancient.vcf.gz.csi",
        vcfstats="popgen/phased.all.info08.name.noancient.stats"
    resources:
        mem_mb=30000,
        runtime=4*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads 4 -S {input.wolves} --force-samples {input.vcf} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """

# ---- FILE FILTER: remove dups ----
rule remove_dups:
    input: 
      vcf="popgen/phased.all.info08.name.noancient.vcf.gz",
      rmfile="king/removedups.txt"
    output: 
      vcf="popgen/phased.all.info08.name.nodups.vcf.gz",
      vcfidx="popgen/phased.all.info08.name.nodups.vcf.gz.csi",
      vcfstats="popgen/phased.all.info08.name.nodups.stats"    
    resources:
      mem_mb=30000,
      runtime=3*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads 4 -S ^{input.rmfile} --force-samples {input.vcf} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """

# ---- FILE FILTER: remove related ----
rule remove_related:
    input: 
      vcf="popgen/phased.all.info08.name.noancient.vcf.gz",
      rmfile="king/removerelated.txt"
    output: 
      vcf="popgen/phased.all.info08.name.norelated.vcf.gz",
      vcfidx="popgen/phased.all.info08.name.norelated.vcf.gz.csi",
      vcfstats="popgen/phased.all.info08.name.norelated.stats"
    params: 
      rmfile="king/removerelated.txt"
    resources:
      mem_mb=30000,
      runtime=3*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads 4 -S ^{input.rmfile} --force-samples {input.vcf} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """

# ---- FILE FILTER: apply MAF filter ----
rule maf_filter:
    input:
      "popgen/phased.all.info08.name.{subset}.vcf.gz"
    output: 
      vcf="popgen/phased.all.info08.name.{subset}.maf01.vcf.gz",
      vcfidx="popgen/phased.all.info08.name.{subset}.maf01.vcf.gz.csi",
      vcfstats="popgen/phased.all.info08.name.{subset}.maf01.stats"
    resources:
      mem_mb=30000,
      runtime=3*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads 4 -i 'MAF>=0.01' {input} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """

# ---- FILE FILTER: European subset (norelated / nodups) ----
rule euro_subset:
    input:
        vcf="popgen/phased.all.info08.name.{subset}.vcf.gz",
        vcfidx="popgen/phased.all.info08.name.{subset}.vcf.gz.csi",
        european=config["eu_wolves_list"]
    output: 
        vcf="popgen/phased.european.info08.name.{subset}.maf01.vcf.gz",
        vcfstats="popgen/phased.european.info08.name.{subset}.maf01.stats"
    resources:
        mem_mb=30000,
        runtime=3*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads 4 -i 'MAF>=0.01' -S {input.european} --force-samples {input.vcf} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """
        
# ---- SMARTPCA: european wolves ----
rule smartpca_subset_european:
    input:
        vcf="popgen/phased.european.info08.name.nodups.maf01.vcf.gz",
        convertf="scripts/convertf_wolves_european.par",
        rscript="scripts/smartpca_related.r",
        smartpca="scripts/smartpca_wolves_european.par",
        euwolves=config["eu_wolves_only_list"]
    output:
        "smartpca/phased.european.info08.name.nodups.maf01.evec",
        "smartpca/phased.european.info08.name.nodups.maf01.eval"
    resources:
        mem_mb = 150000,
        runtime = 10*60
    shell:
        """
        plink --vcf {input.vcf} --dog --double-id --allow-extra-chr \
              --keep-fam {input.euwolves} --recode \
              --out smartpca/phased.european.info08.name.nodups.maf01

        convertf -p {input.convertf}

        Rscript {input.rscript} \
            smartpca/phased.european.info08.name.nodups.maf01.indiv \
            king/removerelated.txt \
            smartpca/phased.european.info08.name.nodups.maf01.norelated.indiv \
            --vanilla

        smartpca -p {input.smartpca}
        """
        
# ---- ADMIXTURE: european wolves ----    
rule admixture_preprocessing_european:
    input:
        vcf="popgen/phased.european.info08.name.norelated.maf01.vcf.gz",
        euwolves=config["eu_wolves_only_list"]
    output:
        vcf="popgen/phased.european.info08.name.norelated.maf01.wolves.vcf.gz",
        vcfstats="popgen/phased.european.info08.name.norelated.maf01.wolves.stats",
        ped="admixture/phased.european.info08.name.norelated.maf01.prune.ped",
        map="admixture/phased.european.info08.name.norelated.maf01.prune.map"
    resources:
        mem_mb = 20000,
        runtime = 3*60
    shell:
        """
        mkdir -p admixture
        # MAF filter and subset for European wolves only
        bcftools view -Oz -i 'MAF>=0.01' -S {input.euwolves} --force-samples {input.vcf} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}

        # Convert to PLINK format
        plink --vcf {output.vcf} \
              --dog --allow-extra-chr --double-id --make-bed \
              --out admixture/phased.european.info08.name.norelated.maf01

        # Fix BIM file naming
        awk '{{ $2 = "chr"$1"_"$4; print }}' \
            admixture/phased.european.info08.name.norelated.maf01.bim > \
            admixture/phased.european.info08.name.norelated.maf01.tmp.bim
        mv admixture/phased.european.info08.name.norelated.maf01.tmp.bim \
           admixture/phased.european.info08.name.norelated.maf01.bim

        # LD prune
        plink --bfile admixture/phased.european.info08.name.norelated.maf01 \
              --indep-pairwise 500 50 0.1 \
              --dog --allow-extra-chr --keep-allele-order \
              --out admixture/phased.european.info08.name.norelated.maf01.prune

        # Recode pruned dataset
        plink --bfile admixture/phased.european.info08.name.norelated.maf01 \
              --extract admixture/phased.european.info08.name.norelated.maf01.prune.prune.in \
              --recode 12 --dog \
              --out admixture/phased.european.info08.name.norelated.maf01.prune
        """

#rule admixture_run_european:
#    input:
#        ped="admixture/phased.european.info08.name.norelated.maf01.prune.ped",
#        map="admixture/phased.european.info08.name.norelated.maf01.prune.map"
#    output:
#        log="admixture/phased.european.info08.name.norelated.maf01.prune_{k}.log",
#        Qfile="admixture/phased.european.info08.name.norelated.maf01.prune_{k}.Q",
#        Pfile="admixture/phased.european.info08.name.norelated.maf01.prune_{k}.P",
#    threads: 4
#    resources:
#        mem_mb = 50000,
#        runtime = 8 * 60
#    shell:
#        """
#        admixture --cv {input.ped} {wildcards.k} -j {threads} \
#            | tee {output.log}
#        mv phased.european.info08.name.norelated.maf01.prune_{k}.Q {output.Qfile}
#        mv phased.european.info08.name.norelated.maf01.prune_{k}.P {output.Pfile}
#        """
        
# ---- Haplonet: admixture european wolves ----
rule haplonet_fileprep:
    input:
        vcf="popgen/phased.european.info08.name.norelated.maf01.wolves.vcf.gz",
    output:
        vcf="haplonet/europeanall/phased.european.info08.name.norelated.maf01.wolves.chr{chr}.vcf.gz"
    resources:
        mem_mb = 150000,
        runtime = 10 * 60
    shell:
        """
        # Split VCF by chromosome
        bcftools view -Oz -r chr{wildcards.chr} {input.vcf} \
            -o {output.vcf}
        bcftools index -f {output.vcf}
        """

rule haplonet_npy:
    input:
        vcf="haplonet/europeanall/phased.european.info08.name.norelated.maf01.wolves.chr{chr}.vcf.gz"
    output:
        npy="haplonet/europeanall/haplonet.chr{chr}.loglike.npy"
    resources:
        mem_mb = 150000,
        runtime = 10 * 60
    conda:
      "haplonet"
    shell:
        """    
        # Train haplonet model for this chromosome
        haplonet train \
            --vcf {input.vcf} --out haplonet/europeanall/haplonet.chr{wildcards.chr}
        """
        
rule haplonet_filelist:
    input:
        npys=expand("haplonet/europeanall/haplonet.chr{chr}.loglike.npy", chr=config["chromosomes"])
    output:
        "haplonet/europeanall/haplonet.filelist"
    run:
        with open(output[0], "w") as f:
            for npy in input.npys:
                f.write(npy + "\n")
                

rule haplonet_admix:
    input:
        filelist = "haplonet/europeanall/haplonet.filelist"
    output:
        q="haplonet/europeanall/haplonet.admixture.k{K}.q",
        npy="haplonet/europeanall/haplonet.admixture.k{K}.f.npy"
    threads: 4
    resources:
        mem_mb = 150000,
        runtime = 15 * 60
    conda:
        "haplonet"
    shell:
        """
        # Run Haplonet admixture
        haplonet admix \
            --filelist {input.filelist} \
            --K {wildcards.K} \
            --threads {threads} \
            --out haplonet/europeanall/haplonet.admixture.k{wildcards.K}
        """
        
rule haplonet_header:
    input:
        vcf = "haplonet/europeanall/phased.european.info08.name.norelated.maf01.wolves.chr38.vcf.gz"
    output:
        list="haplonet/europeanall/names.list"
    resources:
        mem_mb = 20000,
        runtime = 120
    shell:
        """
        bcftools query -l {input.vcf} > {output.list}
        """
                    
# ---- ANGSD: PH ancient wolves for Dstats ----
rule ancient_PH:
    input:
        vcf="popgen/phased.all.info08.name.nodups.maf01.vcf.gz",
        ref=config["ref_genome"],
        bamlist=config["ancient_bams"]
    output:
        outvcf="PH/sites/phased.all.info08.name.nodups.maf01.chr{chr}.positions",
        outvcftv="PH/sites/phased.all.info08.name.nodups.maf01.chr{chr}.tv.positions.txt",
        haplo="PH/angsd_genotypes/{chr}_ancient_PH.haplo.gz",
        plink_bed="PH/angsd_genotypes/{chr}_ancient_PH.bed",
        plink_bim="PH/angsd_genotypes/{chr}_ancient_PH.bim",
        plink_fam="PH/angsd_genotypes/{chr}_ancient_PH.fam",
        plink_positions="PH/sites/plink.{chr}.tv.positions.txt"
    resources:
        mem_mb=20000,
        runtime=4*60
    shell:
        """
        set -euo pipefail
        # extract site list
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' -r chr{wildcards.chr} {input.vcf} -o {output.outvcf}
        # select transversions only
        awk -F '\\t' '(($4 == "A" && $3 == "C") || ($4 == "C" && $3 == "A") || \
        ($4 == "C" && $3 == "G") || ($4 == "G" && $3 == "C") || \
        ($4 == "T" && $3 == "G") || ($4 == "A" && $3 == "T") || \
        ($4 == "T" && $3 == "G") || ($4 == "G" && $3 == "T") || \
        ($4 == "T" && $3 == "A") || $0 ~ /^#/)' {output.outvcf} > {output.outvcftv}
        # index reference
        samtools faidx {input.ref}
        # index sites for ANGSD
        angsd sites index {output.outvcftv}
        # run ANGSD
        angsd -bam {input.bamlist} -checkBamHeaders 0 \
            -dohaplocall 1 -doCounts 1 -doMajorMinor 3 -ref {input.ref} \
            -remove_bads 1 -C 50 -minMapQ 25 -minQ 30 -uniqueOnly 1 -baq 1 \
            -r chr{wildcards.chr} -sites {output.outvcftv} \
            -out PH/angsd_genotypes/{wildcards.chr}_ancient_PH
        # convert to PLINK format
        haploToPlink {output.haplo} PH/angsd_genotypes/{wildcards.chr}_ancient_PH
        # clean zeros in TPED file
        sed -i 's/N/0/g' PH/angsd_genotypes/{wildcards.chr}_ancient_PH.tped
        # make PLINK binary files
        plink --tfile PH/angsd_genotypes/{wildcards.chr}_ancient_PH \
            --allow-extra-chr --dog --make-bed --double-id \
            --out PH/angsd_genotypes/{wildcards.chr}_ancient_PH --keep-allele-order
        awk '{{ print $1 "_" $2 }}' {output.outvcftv} > {output.plink_positions}       
        """
        
rule ancient_PH_combine:
    input:
        bed=expand("PH/angsd_genotypes/{chr}_ancient_PH.bed", chr=config["chromosomes"]),
        bim=expand("PH/angsd_genotypes/{chr}_ancient_PH.bim", chr=config["chromosomes"]),
        fam=expand("PH/angsd_genotypes/{chr}_ancient_PH.fam", chr=config["chromosomes"]),
        plink_positions=expand("PH/sites/plink.{chr}.tv.positions.txt", chr=config["chromosomes"]),
        modern_bed="Dstats/phased.all.info08.name.nodups.maf01.bed",
        modern_bim="Dstats/phased.all.info08.name.nodups.maf01.bim",
        modern_fam="Dstats/phased.all.info08.name.nodups.maf01.fam",
        bamnames="data/ancientwolves_bamnames.txt",
        convertfancient="scripts/convertf_ancient.par",
        convertfmodern="scripts/convertf_modern.par",
        mergeit="scripts/mergeit_ancient_modern.par"
    output:
        combined_ped="PH/angsd_genotypes/combined.ped",
        combined_map="PH/angsd_genotypes/combined.map",
        fam="PH/phased.all.info08.name.nodups.maf01.tv.fam",
        bed="PH/phased.all.info08.name.nodups.maf01.tv.bed",
        bim="PH/phased.all.info08.name.nodups.maf01.tv.bim",
        merged_ind="PH/ancient_modern_merged.ind",
        tv_sites="PH/sites/plink.tv.positions.txt"
    resources:
        mem_mb=60000,
        runtime=4*60
    shell:
        """
        # Combine all chromosome position files
        rm -f {output.tv_sites}
        cat {input.plink_positions} >> {output.tv_sites}

        # Create list of ancient files
        rm -f PH/angsd_genotypes/filelist
        for f in {input.bed}; do
            echo "${{f%.bed}}" >> PH/angsd_genotypes/filelist
        done

        # Merge ancient chromosomes
        plink --merge-list PH/angsd_genotypes/filelist \
              --allow-extra-chr --dog --recode \
              --out PH/angsd_genotypes/combined
        
        # Replace names in ped file
        awk 'NR==FNR {{a[NR]=$1; next}} {{$1=a[FNR]; $2=a[FNR]; print}}' \
        {input.bamnames} PH/angsd_genotypes/combined.ped > PH/angsd_genotypes/combined.tmp.ped
        mv PH/angsd_genotypes/combined.tmp.ped PH/angsd_genotypes/combined.ped
              
        awk '{{ $2 = "chr" $1 "_" $4; print }}' Dstats/phased.all.info08.name.nodups.maf01.bim > Dstats/phased.all.info08.name.nodups.maf01.snp.bim
        mv Dstats/phased.all.info08.name.nodups.maf01.snp.bim Dstats/phased.all.info08.name.nodups.maf01.bim
        
        # Extract TV sites from modern
        plink --bfile Dstats/phased.all.info08.name.nodups.maf01 \
              --extract {output.tv_sites} \
              --dog --allow-extra-chr --recode \
              --out PH/phased.all.info08.name.nodups.maf01.tv

        plink --bfile Dstats/phased.all.info08.name.nodups.maf01 \
              --extract {output.tv_sites} \
              --dog --allow-extra-chr --make-bed \
              --out PH/phased.all.info08.name.nodups.maf01.tv

        # Convert to eigenstrat
        convertf -p {input.convertfmodern}
        convertf -p {input.convertfancient}

        # Merge ancient + modern
        mergeit -p {input.mergeit}
        """

        
# ---- ADMIXTOOLS: Dstats and f4ratio ----  
rule dstat_plink_convert:
    input:
        vcf="popgen/phased.all.info08.name.nodups.maf01.vcf.gz"
    output:
        bed="Dstats/phased.all.info08.name.nodups.maf01.bed",
        bim="Dstats/phased.all.info08.name.nodups.maf01.bim",
        fam="Dstats/phased.all.info08.name.nodups.maf01.fam"
    resources:
        mem_mb=100000,
        runtime=5*60
    shell:
        """
        # Convert to PLINK format
        plink --vcf {input.vcf} \
              --dog --double-id --allow-extra-chr \
              --make-bed \
              --out Dstats/phased.all.info08.name.nodups.maf01

        # Fix BIM file naming
        awk '{{ $2 = "chr" $1 "_" $4; print }}' \
            Dstats/phased.all.info08.name.nodups.maf01.map \
            > Dstats/phased.all.info08.name.nodups.maf01.snp.map

        mv Dstats/phased.all.info08.name.nodups.maf01.snp.map \
           Dstats/phased.all.info08.name.nodups.maf01.map
        """
        
rule dstat_dog:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/admixtools_Dstat_dog.r"
    output:
          nwib="Dstats/dog/Dstat_dog_NWIB.list",
          scan="Dstats/dog/Dstat_dog_SCAN.list",
          dstats1="Dstats/dog/Dstat_dog_NWIB.rds",
          dstats2="Dstats/dog/Dstat_dog_SCAN.rds"
    conda:
        "Renv"
    resources:
        mem_mb=50000,
        runtime=6*60
    shell:
        """
        # MW122 (NWIB)
        awk '{{print "AndeanFox\tBasenjiDog\tMW122\t" $1}}' {input.indiv} > {output.nwib}
        awk '{{print "AndeanFox\tTarabaDog\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tAfghanDog\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tEG44\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tHebeiDog\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tGansu2Dog\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tAlaskanMalamuteDog\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tTasiilaq_51603\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tNovembre_Dingo\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tPG84\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tGShepDog\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tFinnishLapphundDog\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        # MW145 (SCAN)
        awk '{{print "AndeanFox\tBasenjiDog\tMW145\t" $1}}' {input.indiv} > {output.scan}
        awk '{{print "AndeanFox\tTarabaDog\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tAfghanDog\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tEG44\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tHebeiDog\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tGansu2Dog\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tAlaskanMalamuteDog\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tTasiilaq_51603\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tNovembre_Dingo\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tPG84\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tGShepDog\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tFinnishLapphundDog\tMW145\t" $1}}' {input.indiv} >> {output.scan}

        Rscript {input.scripts} --vanilla
        """
rule dstat_wild:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/admixtools_Dstat_wild.r"
    output:
        nwib = "Dstats/wild/Dstat_wild_NWIB.list",
        scan = "Dstats/wild/Dstat_wild_SCAN.list",
        dstats1 = "Dstats/wild/Dstat_wild_NWIB.rds",
        dstats2 = "Dstats/wild/Dstat_wild_SCAN.rds"
    conda:
        "Renv"
    resources:
        mem_mb=50000,
        runtime=6*60
    shell:
        """
        # MW122 (NWIB)
        awk '{{print "AndeanFox\tAlgerian_grey_wolf\tMW122\t" $1}}' {input.indiv} > {output.nwib}
        awk '{{print "AndeanFox\tBerlinZoo\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tC_adustus\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tC_mesomelas\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tCaliforniaCoyote\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tEthiopian_grey_wolf\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tKenyaJackal\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tKruger\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tSyrian_jackal\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tZimbabwe1\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        # MW145 (SCAN)
        awk '{{print "AndeanFox\tAlgerian_grey_wolf\tMW145\t" $1}}' {input.indiv} > {output.scan}
        awk '{{print "AndeanFox\tBerlinZoo\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tC_adustus\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tC_mesomelas\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tCaliforniaCoyote\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tEthiopian_grey_wolf\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tKenyaJackal\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tKruger\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tSyrian_jackal\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tZimbabwe1\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        
        Rscript {input.scripts} --vanilla
        """
rule dstat_ancient:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/admixtools_Dstat_ancient.r"
    output:
        nwib = "Dstats/ancient/Dstat_ancient_NWIB.list",
        scan = "Dstats/ancient/Dstat_ancient_SCAN.list",
        dstats1="Dstats/ancient/Dstat_ancient_SCAN.rds",
        dstats2="Dstats/ancient/Dstat_ancient_NWIB.rds"
    conda:
        "Renv"
    resources:
        mem_mb=50000,
        runtime=6*60
    shell:
        """
        # MW122 (NWIB)
        awk '{{print "AndeanFox\tIN18-005\tMW122\t" $1}}' {input.indiv} > {output.nwib}
        awk '{{print "AndeanFox\tAL2657\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tCGG33\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tJK2181\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tTumat2\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tJK2179\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tTU839\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tTU840\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tAL2370\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tDS-CANIS-HMNH-011\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tAL3185\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tAL2350\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tPON012\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tCANIS-HMNH-007\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        # MW145 (SCAN)
        awk '{{print "AndeanFox\tIN18-005\tMW145\t" $1}}' {input.indiv} > {output.scan}
        awk '{{print "AndeanFox\tAL2657\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tCGG33\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tJK2181\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tTumat2\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tJK2179\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tTU839\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tTU840\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tAL2370\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tDS-CANIS-HMNH-011\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tAL3185\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tAL2350\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tPON012\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tCANIS-HMNH-007\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        Rscript {input.scripts} --vanilla
        """
rule dstat_asian:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/admixtools_Dstat_asian.r"
    output:
        nwib = "Dstats/asian/Dstat_asian_NWIB.list",
        scan = "Dstats/asian/Dstat_asian_SCAN.list",
        dstats1="Dstats/asian/Dstat_asian_SCAN.rds",
        dstats2="Dstats/asian/Dstat_asian_NWIB.rds"
    conda:
        "Renv"
    resources:
        mem_mb=50000,
        runtime=6*60
    shell:
        """
        # MW122 (NWIB)
        awk '{{print "AndeanFox\tMW1063\tMW122\t" $1}}' {input.indiv} > {output.nwib}
        awk '{{print "AndeanFox\tMW522\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        awk '{{print "AndeanFox\tMW167\tMW122\t" $1}}' {input.indiv} >> {output.nwib}
        # MW145 (SCAN)
        awk '{{print "AndeanFox\tMW1063\tMW145\t" $1}}' {input.indiv} > {output.scan}
        awk '{{print "AndeanFox\tMW522\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        awk '{{print "AndeanFox\tMW167\tMW145\t" $1}}' {input.indiv} >> {output.scan}
        Rscript {input.scripts} --vanilla
        """

rule dstat_asianregions:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/admixtools_Dstat_asianregions.r"
    output:
        list="Dstats/asianregions/Dstat_asianregions.list",
        dstats="Dstats/asianregions/Dstat_asianregions.rds"
    conda:
        "Renv"
    resources:
        mem_mb=50000,
        runtime=6*60
    shell:
        """
        awk '{{print "AndeanFox\t" $1 "\tMW1076\tMW497"}}' {input.indiv} > {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN16\tMW497"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW530\tMW497"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW564\tMW497"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW535\tMW497"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN9\tMW497"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tKoreanWolf4\tMW497"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN16\tMW1076"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW530\tMW1076"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW564\tMW1076"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW535\tMW1076"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN9\tMW1076"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tKoreanWolf4\tMW1076"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW530\tChinese_CAN16"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW564\tChinese_CAN16"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW535\tChinese_CAN16"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN9\tChinese_CAN16"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tKoreanWolf4\tChinese_CAN16"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW564\tMW530"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW535\tMW530"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN9\tMW530"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tKoreanWolf4\tMW530"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tMW535\tMW564"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN9\tMW564"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tKoreanWolf4\tMW564"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tChinese_CAN9\tMW535"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tKoreanWolf4\tMW535"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tKoreanWolf4\tChinese_CAN9"}}' {input.indiv} >> {output.list}
        Rscript {input.scripts} --vanilla
        """
        
rule dstat_ancientdog:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/admixtools_Dstat_ancientdog.r"
    output:
        list="Dstats/ancientdog/Dstat_ancientdog.list",
        dstats="Dstats/ancientdog/Dstat_ancientdog.rds"
    conda:
        "Renv"
    resources:
        mem_mb=50000,
        runtime=6*60
    shell:
        """
        awk '{{print "AndeanFox\t" $1 "\tAL2657\tGShepDog"}}' {input.indiv} > {output.list}
        awk '{{print "AndeanFox\t" $1 "\tCGG33\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tJK2181\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tTumat2\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tJK2179\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tTU839\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tTU840\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAL2370\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tDS-CANIS-HMNH-011\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAL3185\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAL2350\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tPON012\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tCANIS-HMNH-007\tGShepDog"}}' {input.indiv} >> {output.list}
        Rscript {input.scripts} --vanilla
        """
             
rule dstat_dogdog:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/admixtools_Dstat_dogdog.r"
    output:
        list="Dstats/dogdog/Dstat_dogdog.list",
        dstats="Dstats/dogdog/Dstat_dogdog.rds"
    conda:
        "Renv"
    resources:
        mem_mb=50000,
        runtime=6*60
    shell:
        """
        awk '{{print "AndeanFox\t" $1 "\tBasenjiDog\tGShepDog"}}' {input.indiv} > {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAlaskanHusky_SY001\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tNovembre_Dingo\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tHebeiDog\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAfghanDog\tGShepDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tNovembre_Dingo\tBasenjiDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAlaskanHusky_SY001\tBasenjiDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tHebeiDog\tBasenjiDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAfghanDog\tBasenjiDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tHebeiDog\tNovembre_Dingo"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAlaskanHusky_SY001\tNovembre_Dingo"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAfghanDog\tNovembre_Dingo"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAfghanDog\tHebeiDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAlaskanHusky_SY001\tHebeiDog"}}' {input.indiv} >> {output.list}
        awk '{{print "AndeanFox\t" $1 "\tAfghanDog\tAlaskanHusky_SY001"}}' {input.indiv} >> {output.list}
        Rscript {input.scripts} --vanilla
        """

rule f4_ratio:
    input:
        indiv = "PH/phased.all.info08.name.nodups.maf01.tv.fam",
        scripts = "scripts/f4ratio_1.r"
    output:
        "qpadm/f4ratio_all.RData"
    conda:
        "Renv"
    resources:
        mem_mb=60000,
        runtime=24*60
    shell:
        """
        Rscript {input.scripts} --vanilla
        """

        
# ---- ORIENTAGRAPH: topology  ----    
rule orientagraph_preprocessing_all:
    input:
        ids="scripts/europeanwolves_Afox.txt",
        bed="popgen/phased.all.info08.name.bed",
        bim="popgen/phased.all.info08.name.bim",
        fam="popgen/phased.all.info08.name.fam",
        snps="admixture/phased.european.info08.name.norelated.maf01.prune.prune.in",
        renamepops="scripts/populations.r"
    params:
        plink="popgen/phased.all.info08.name"
    output:
        strat="orientagraph/phased.european.fox.info08.input.frq.strat.gz",
        bim="orientagraph/phased.european.fox.info08.name.prune.bim"
    resources:
        mem_mb = 20000,
        runtime = 2*60
    conda:
        "Renv" 
    shell:
        """
        #rename the snps in the plink file
        awk '{{ $2 = "chr" $1 "_" $4; print }}' {params.plink}.bim > popgen/phased.all.info08.name.snp.bim
        mv popgen/phased.all.info08.name.snp.bim {params.plink}.bim 
        #filter for European wolves
        plink --bfile {params.plink} --keep-fam {input.ids} --extract {input.snps} --make-bed --dog \
        --out orientagraph/phased.european.fox.info08.name.prune
        #rename the populations
        Rscript {input.renamepops}  --vanilla
        #convert to the treemix format
        plink --bfile orientagraph/phased.european.fox.info08.name.prune --family --dog --freq --out orientagraph/phased.european.fox.info08.input
        gzip orientagraph/phased.european.fox.info08.input.frq.strat
        """

rule orientagraph_preprocessing_all2:
    input:
        strat="orientagraph/phased.european.fox.info08.input.frq.strat.gz",
        py="scripts/plink2treemix.py"
    output:
        ogfile="orientagraph/phased.european.fox.info08.input.frq.gz"
    resources:
        mem_mb = 20000,
        runtime = 1*60
    shell:
        """
        module purge
        module load python/2.7.17
        module load plink/1.9.0
        module load openblas/0.3.24
        python "{input.py}" "{input.strat}" "{output.ogfile}" || echo "Script exited with $?"
        """        

rule run_orientagraph_all:
    input:
        freq="orientagraph/phased.european.fox.info08.input.frq.gz" 
    output:
        "orientagraph/og_rep/allpop_m{m}_rep{r}.treeout.gz",
        log="orientagraph/og_log/log_allpop_m{m}_rep{r}.log"
    resources:
        mem_mb=80000,
        runtime=10*60
    conda:
        "orientagraph"
    shell:
        """
        set -euo pipefail
        # Run OrientAGraph
        orientagraph -i {input.freq} -m {wildcards.m} \
        -o orientagraph/og_rep/allpop_m{wildcards.m}_rep{wildcards.r} \
        -k 2000 -root "AndeanFox" -mlno > {output.log} 2>&1
        # Fail if treeout.gz was not created
        if [ ! -f {output[0]} ]; then
            echo "Error: {output[0]} not created" >&2
            exit 1
        fi
        """

rule orientagraph_preprocessing_noSWASIA:
    input:
        bim="orientagraph/phased.european.fox.info08.name.prune.bim",
        py="scripts/plink2treemix.py"
    params:
      plink="orientagraph/phased.european.fox.info08.name.prune"
    output:
        ogfile="orientagraph/phased.european.fox.info08.name.prune.noSWASIA.input.frq.gz"
    resources:
        mem_mb = 20000,
        runtime = 2*60
    shell:
        """
        module purge
        module load python/2.7.17
        module load plink/1.9.0
        module load openblas/0.3.24
        echo "SWASIA" > orientagraph/SWASIA.txt
        plink --bfile {params.plink} --remove-fam orientagraph/SWASIA.txt --make-bed --dog \
        --out orientagraph/phased.european.fox.info08.name.prune.noSWASIA    

        #convert to the treemix format
        plink --bfile orientagraph/phased.european.fox.info08.name.prune.noSWASIA --family --dog --freq \
        --out orientagraph/phased.european.fox.info08.name.prune.noSWASIA.input
        gzip orientagraph/phased.european.fox.info08.name.prune.noSWASIA.input.frq.strat
        python {input.py} orientagraph/phased.european.fox.info08.name.prune.noSWASIA.input.frq.strat.gz  {output.ogfile}
        """

rule run_orientagraph_noSWASIA:
    input:
        freq="orientagraph/phased.european.fox.info08.name.prune.noSWASIA.input.frq.gz" 
    output:
        "orientagraph/og_rep2/allpop_m{m}_rep{r}.treeout.gz",
        log="orientagraph/og_log2/log_allpop_m{m}_rep{r}.log"
    resources:
        mem_mb=120000,
        runtime=15*60
    conda:
        "orientagraph"
    shell:
        """
        set -euo pipefail
        # Run OrientAGraph
        orientagraph -i {input.freq} -m {wildcards.m} \
            -o orientagraph/og_rep2/allpop_m{wildcards.m}_rep{wildcards.r} \
            -k 2000 -root "AndeanFox" -mlno \
            > {output.log} 2>&1

        # Fail if treeout.gz was not created
        if [ ! -f {output[0]} ]; then
            echo "Error: {output[0]} not created" >&2
            exit 1
        fi
        """

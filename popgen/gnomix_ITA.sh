#!/bin/bash
#SBATCH --job-name=gnomixITA
#SBATCH --output=gnomixITA.out
#SBATCH --error=gnomixITA.err
#SBATCH --ntasks=1
#SBATCH --array=1-38              
#SBATCH --cpus-per-task=1         
#SBATCH --mem-per-cpu=20G         
#SBATCH --time=00:40:00           

###############################################################################
# Author: Evelyn Todd
# Date: 26/08/2025
#
# Local ancestry inference using Gnomix for ITA wolf hybrids.
###############################################################################

mkdir -p gnomix/ITA_chr${SLURM_ARRAY_TASK_ID}
mkdir -p gnomix/positions
mkdir -p gnomix/map
mkdir -p gnomix/ITA
mkdir -p gnomix/ITAall

source ~/.bashrc

module load gsl/2.5 perl/5.38.0 bcftools/1.20

conda activate gnomix

###############################################################################
# Create chromosome positions from the pruned ADMIXTURE dataset
###############################################################################

awk '{print "chr"$1 "\t" $4}' admixture/combined.europeanwolves.norelated.maf01.prune.map | grep "chr${SLURM_ARRAY_TASK_ID}" > gnomix/positions/positions.chr${SLURM_ARRAY_TASK_ID}.txt

###############################################################################
# Create target and reference VCFs
###############################################################################

# Extract hybrid individuals
bcftools view -S data/MOSAICpops/ITA_hybrid.txt -R gnomix/positions/positions.chr${SLURM_ARRAY_TASK_ID}.txt popgen/combined.nodups.vcf.gz -o gnomix/ITA_chr${SLURM_ARRAY_TASK_ID}/target_chr${SLURM_ARRAY_TASK_ID}.vcf

# Create reference sample list
cat data/MOSAICpops/Dog.txt > data/MOSAICpops/ITADog.txt
cat data/MOSAICpops/ITA.txt >> data/MOSAICpops/ITADog.txt

# Extract reference individuals
bcftools view -S data/MOSAICpops/ITADog.txt -R gnomix/positions/positions.chr${SLURM_ARRAY_TASK_ID}.txt popgen/combined.nodups.vcf.gz -o gnomix/ITA_chr${SLURM_ARRAY_TASK_ID}/ref_chr${SLURM_ARRAY_TASK_ID}.vcf

###############################################################################
# Prepare Gnomix input files
###############################################################################

# Metadata file containing reference population labels. Edit manually
# cp /bin/gnomix/config.yaml gnomix/ITA/config.yaml

# Create genetic map using canfam recombination map
awk 'BEGIN{OFS="\t"; print "chm","pos","pos_cm"}
NR==1 {next}
$1 ~ /^Chromosome/ {next}
{print "chr"$1, $2, $4}' map/chr${SLURM_ARRAY_TASK_ID}_map.txt > gnomix/map/chr${SLURM_ARRAY_TASK_ID}.map

# Format sample map
awk 'BEGIN{OFS="\t"} {print $1,$2}' gnomix/ITA/reftest.smap > gnomix/ITA/refITA.smap

###############################################################################
# Run Gnomix
###############################################################################

python3 gnomix.py \
    gnomix/ITA_chr${SLURM_ARRAY_TASK_ID}/target_chr${SLURM_ARRAY_TASK_ID}.vcf \
    gnomix/ITA_chr${SLURM_ARRAY_TASK_ID} \
    chr${SLURM_ARRAY_TASK_ID} \
    False \
    gnomix/map/chr${SLURM_ARRAY_TASK_ID}.map \
    gnomix/ITA_chr${SLURM_ARRAY_TASK_ID}/ref_chr${SLURM_ARRAY_TASK_ID}.vcf \
    gnomix/ITA/refITA.smap \
    gnomix/ITA/config.yaml

###############################################################################
# Save output
###############################################################################

cp gnomix/ITA_chr${SLURM_ARRAY_TASK_ID}/query_results.lai gnomix/ITAall/query_results_chr${SLURM_ARRAY_TASK_ID}.lai
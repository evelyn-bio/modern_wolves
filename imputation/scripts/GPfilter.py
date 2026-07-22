#!/usr/bin/env python3
"""
Filter variants by GP (genotype posterior) quality across samples.
Outputs BED file of variants where =X% of samples have max(GP) = threshold.
"""
import sys
from cyvcf2 import VCF

# Parse arguments
if len(sys.argv) != 5:
    print("Usage: python filter_by_gp.py <input.vcf.gz> <output.bed> <gp_threshold> <sample_threshold>")
    print("Example: python filter_by_gp.py input.bcf output.bed 0.9 0.9")
    sys.exit(1)

vcf_in = sys.argv[1]
out_file = sys.argv[2]
gp_threshold = float(sys.argv[3])      # 0.9,0.95 etc
sample_threshold = float(sys.argv[4])  # 0.9 (= 90% of samples)

# Open VCF/BCF
vcf = VCF(vcf_in)
n_samples = len(vcf.samples)

print(f"Processing {vcf_in}")
print(f"Number of samples: {n_samples}")
print(f"GP threshold: {gp_threshold}")
print(f"Sample threshold: {sample_threshold} ({sample_threshold*100}% of samples)")
print(f"Minimum passing samples: {int(n_samples * sample_threshold)}")

# Output file
out = open(out_file, "w")

# Counters
total_variants = 0
passed_variants = 0

# Process variants
for v in vcf:
    total_variants += 1
    
    # Progress indicator
    if total_variants % 10000 == 0:
        print(f"Processed {total_variants} variants, kept {passed_variants}", file=sys.stderr)
    
    # Get GP field
    gp = v.format("GP")
    
    # Skip if GP field is missing
    if gp is None:
        continue
    
    # Count samples passing GP threshold
    passed = 0
    for sample_gp in gp:
        # Handle None or malformed values
        if sample_gp is None:
            continue
        
        # Check if max GP >= threshold
        # sample_gp is an array of probabilities for each genotype
        if max(sample_gp) >= gp_threshold:
            passed += 1
    
    # Calculate fraction
    frac = passed / n_samples
    
    # Write to BED if passes sample threshold
    if frac >= sample_threshold:
        # BED format: chrom, start (0-based), end (1-based)
        # VCF POS is 1-based, so BED start = POS-1, end = POS
        out.write(f"{v.CHROM}\t{v.POS}\n")
        passed_variants += 1

out.close()

# Summary
print(f"\nDone!")
print(f"Total variants: {total_variants}")
print(f"Passed variants: {passed_variants}")
print(f"Pass rate: {passed_variants/total_variants*100:.2f}%")
print(f"Output written to: {out_file}")
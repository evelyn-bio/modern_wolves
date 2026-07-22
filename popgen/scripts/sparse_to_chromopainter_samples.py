#!/usr/bin/env python3

import gzip
import re
import sys

if len(sys.argv) != 3:
    print("Usage:")
    print("python sparse_to_chromopainter_samples.py input.samples.txt.gz output.samples.out.gz")
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

# ----------------------------------------------------------------------
# First pass: count donor haplotypes (= number of painting samples)
# ----------------------------------------------------------------------

nhaps = 0

with gzip.open(infile, "rt") as f:
    for line in f:
        line = line.strip()
        if re.match(r".+_[01]$", line):
            nhaps += 1

# ----------------------------------------------------------------------
# Second pass: rewrite into ChromoPainter format
# ----------------------------------------------------------------------

nind = nhaps // 2

header = (
    f"EM_iter = 0 (N_e = 0 / copy_prop = 0 / mutation = 0 / "
    f"mutationGLOBAL = 0), nsamples = {nind}\n"
)

with gzip.open(infile, "rt") as fin, gzip.open(outfile, "wt") as fout:

    fout.write(header)

    for line in fin:

        line = line.rstrip()

        m = re.match(r"(.+)_([01])$", line)

        if m:
            sample = m.group(1)
            hap = int(m.group(2)) + 1

            fout.write(f"HAP {hap} {sample}\n")
        else:
            fout.write(line + "\n")

print(f"Wrote {outfile}")
print(f"Detected {nhaps} donor haplotypes.")
#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

if len(sys.argv) != 4:
    print("Usage:")
    print("python mosaic_map.py <genetic_map> <positions_file> <output.rec>")
    sys.exit(1)

mapfile = sys.argv[1]
posfile = sys.argv[2]
outfile = sys.argv[3]

# Read genetic map
gmap = pd.read_csv(mapfile, sep=r"\s+")

# Read SNP positions (chr pos)
pos = pd.read_csv(
    posfile,
    sep=r"\s+",
    header=None,
    names=["chr", "pos"]
)

# Interpolate cumulative genetic map
map_cm = np.interp(
    pos["pos"],
    gmap["Position(bp)"],
    gmap["Map(cM)"]
)

# Write MOSAIC rec file
with open(outfile, "w") as out:
    out.write(f":sites:{len(pos)}\n")
    out.write(" ".join(map(str, pos["pos"])) + "\n")
    out.write(" ".join(f"{x:.10f}" for x in map_cm) + "\n")

print(f"Wrote {outfile}")
print(f"Interpolated {len(pos)} SNPs.")
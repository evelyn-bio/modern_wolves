#!/usr/bin/env python3

import sys

with open(sys.argv[1]) as f:
    rows = [line.rstrip() for line in f]

with open(sys.argv[2], "w") as out:
    for col in zip(*rows):
        out.write("".join(col) + "\n")
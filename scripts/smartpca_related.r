#!/usr/bin/env Rscript

# Remove duplicate IDs from the Dstat input file

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: script.R <indiv_file> <related_file> <output_file>")
}

indiv_file   <- args[1]
related_file <- args[2]
output_file  <- args[3]

# Read data
indiv <- read.table(indiv_file, header = FALSE)
related <- read.table(related_file, header = FALSE)

# Mark related/unrelated
indiv$V3 <- ifelse(indiv$V1 %in% related$V1, "related", "unrelated")

# Optional: print summary info
str(indiv)
str(related)
table(indiv$V3)

# Write output
write.table(indiv, output_file,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t")

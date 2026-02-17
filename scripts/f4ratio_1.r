##Author: Evelyn Todd
##Date: 27/09/2025
##Description: qpadm for dog introgression
setwd("/projects/psg/people/pkb156/MW")
library(admixtools)
library(tidyverse)

##################################################
##read in the names of all individuals to test-----------------
df<-read.table("PH/phased.all.info08.name.nodups.maf01.tv.fam", header=F)
str(df) 

##################################################
## Build models for qpf4ratio --------------------
# qpf4ratio expects a 5-column matrix: p1,p2,p3,p4,p5
# Each row is one ratio to compute.
models <- expand_grid(
        pop1 = "FinnishLapphundDog",
        pop2 = "AndeanFox",
        pop3 = unique(df$V1),
        pop4 = "MW122",
        pop5 = "GShepDog")

head(models)


prefix = 'PH/phased.all.info08.name.nodups.maf01.tv'

batch_size <- 100  # tune depending on memory
n_models <- nrow(models)
n_batches <- ceiling(n_models / batch_size)

res_list <- vector("list", n_batches)

for (i in seq_len(n_batches)) {
  start <- (i - 1) * batch_size + 1
  end <- min(i * batch_size, n_models)
  cat("Running batch", i, "of", n_batches, "(", start, "-", end, ")\n")
  
  batch <- models[start:end, ]
  res <- qpf4ratio(data = prefix, pops = as.matrix(batch))
  res_list[[i]] <- res
  
  # Optional: save intermediate results
  save(res, file = sprintf("qpadm/f4ratio_batch_%02d.RData", i))
  
  rm(res); gc()
}

# Combine all
res_f4 <- do.call(rbind, res_list)
save(res_f4, file = "qpadm/f4ratio_all.RData")
#Description: Dstats of ancient admixture into wolves
#Author: Evelyn Todd
#Date: 15/09/25
####################
##set env ----
setwd("/projects/psg/people/pkb156/MW")
library(admixtools)
df1<-read.table("Dstats/ancientdog/Dstat_ancientdog.list", header=F)

Dstats1<-qpdstat("PH/ancient_modern_merged",
        pop1 = unique(df1$V1),
        pop2 = unique(df1$V2),
        pop3 = unique(df1$V3),
        pop4 = unique(df1$V4),
        unique_only=TRUE,
        allsnps = TRUE,
        f4mode=FALSE)
str(Dstats1)

saveRDS(Dstats1, file = "Dstats/ancientdog/Dstat_ancientdog.rds")

#Description: Dstats of dog admixture into wolves
#Author: Evelyn Todd
#Date: 12/09/25
####################
##set env ----
setwd("/projects/psg/people/pkb156/MW")
library(admixtools)
df1<-read.table("Dstats/dog/Dstat_dog_NWIB.list", header=F)

Dstats1<-qpdstat("PH/phased.all.info08.name.nodups.maf01.tv",
        pop1 = unique(df1$V1),
        pop2 = unique(df1$V2),
        pop3 = unique(df1$V3),
        pop4 = unique(df1$V4),
        unique_only=TRUE,
        f4mode=FALSE)
str(Dstats1)

saveRDS(Dstats1, file = "Dstats/dog/Dstat_dog_NWIB.rds")

df2<-read.table("Dstats/dog/Dstat_dog_SCAN.list", header=F)

Dstats2<-qpdstat("PH/phased.all.info08.name.nodups.maf01.tv",
        pop1 = unique(df2$V1),
        pop2 = unique(df2$V2),
        pop3 = unique(df2$V3),
        pop4 = unique(df2$V4),
        unique_only=TRUE,
        f4mode=FALSE)
str(Dstats2)

saveRDS(Dstats2, file = "Dstats/dog/Dstat_dog_SCAN.rds")
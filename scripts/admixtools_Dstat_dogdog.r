#Description: Dstats of ancient admixture into wolves
#Author: Evelyn Todd
#Date: 15/09/25
####################
##set env ----
setwd("/projects/psg/people/pkb156/MW")
library(admixtools)
df1<-read.table("Dstats/dogdog/Dstat_dogdog.list", header=F)

Dstats1<-qpdstat("PH/phased.all.info08.name.nodups.maf01.tv",
        pop1 = unique(df1$V1),
        pop2 = unique(df1$V2),
        pop3 = unique(df1$V3),
        pop4 = unique(df1$V4),
        unique_only=TRUE,
        allsnps = TRUE,
        f4mode=FALSE)
str(Dstats1)

saveRDS(Dstats1, file = "Dstats/dogdog/Dstat_dogdog.rds")

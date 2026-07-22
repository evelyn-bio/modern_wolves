#Description: Dstats of ancient admixture into wolves
#Author: Evelyn Todd
#Date: 15/09/25
####################
##set env ----
setwd("/projects/psg/people/pkb156/MW_sub2")
library(admixtools)
library(data.table)
library(tidyverse)
df1<-read.table("Dstats/ancientmodern/Dstat_ancientmodern.list", header=F)

meta<- fread("data/sampleinformation_corrected.txt", header=T)
str(meta)
meta<-subset(meta, meta$relatedness =="unrelated")
euwolves<-c("CEL", "CAR", "DINBAL","ITA", "NWIB")
df2 <- df1 %>%
  left_join(meta %>% select(ID, newgroup), by = c("V3" = "ID")) %>%
  rename(pop3group = newgroup) %>%
  left_join(meta %>% select(ID, newgroup), by = c("V4" = "ID")) %>%
  rename(pop4group = newgroup) %>% 
  filter(pop4group != "NA") %>% 
  filter(pop3group != "NA") %>% 
  filter(pop3group != pop4group)  %>%
  filter(pop3group %in% c(c(euwolves))) %>%
  filter(pop4group %in% c(c(euwolves)))

Dstats1<-qpdstat("PH/ancient_modern_merged",
        pop1 = unique(df2$V1),
        pop2 = unique(df2$V2),
        pop3 = unique(df2$V3),
        pop4 = unique(df2$V4),
        unique_only=TRUE,
        allsnps = TRUE,
        f4mode=FALSE)
str(Dstats1)

saveRDS(Dstats1, file = "Dstats/ancientmodern/Dstat_ancientmodern.rds")

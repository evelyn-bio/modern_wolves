#Author: Evelyn Todd
#Date: 06.10.25
#Description: add family groups to fam file for orientagraph
setwd("/projects/psg/people/pkb156/MW") 
meta<-read.table("data/allcanids_metadata_ancestrygroups_edited_nonrelated_updated.txt",header=T, sep="\t") 
fam<-read.table("orientagraph/phased.european.fox.info08.name.prune.fam") 
str(meta) 
str(fam) 
merged<-merge(fam, meta, by.x="V1", by.y="ID", all.x=T, sort=F)
merged$newgroup<-ifelse(is.na(merged$newgroup),merged$V1, merged$newgroup) 
final<-merged[, c("newgroup", "V1", "V3", "V4", "V5", "V6")] 
str(final)
write.table(final, "orientagraph/phased.european.fox.info08.name.prune.fam", sep="\t", quote=F, row.names=F, col.names=F)

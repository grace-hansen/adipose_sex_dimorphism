#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
arguments <- parse_args(OptionParser(usage = "plot_exp.R <trait>",option_list=list()),
                        positional_arguments = 1)
opt<-arguments$opt
trait<-arguments$args[1]

setwd(paste("/home/grace/midway/ldsc_seg",sep=""))
pony_colors<-fread("~/papers/TWAS/pony_palette")


dat<-fread(paste("output/",trait,"/GTEx_exp/all_tissues.results",sep=''))
groups<-fread(paste("output/",trait,"/GTEx_exp/tissues_groups",sep=''))
dat<-merge(groups,dat,by="Tissue")
dat<-dat %>% arrange(Group,-log10(Enrichment_p))
dat<-dat %>% arrange(Group,Enrichment)
dat$Tissue<-factor(dat$Tissue,levels=unique(dat$Tissue))

set.seed(29)
pony_palette=slice(pony_colors,sample(1:length(unique(dat$Group))))

G<-ggplot(dat,aes(x=Tissue,y=-log10(Enrichment_p)))+
  geom_bar(aes(fill=Group),position="dodge",stat="identity")+
  theme_minimal()+
  scale_y_continuous("-log10 Enrichment p-value")+
  scale_fill_manual(values=rgb(pony_palette[,1:3]))+
  theme(axis.title=element_text(size=15),
        legend.position = "none")+
  scale_x_discrete(breaks = NULL)
G

pdf(paste("output/",trait,"/GTEx_exp/GTEx_ldscseg_exp.pdf",sep=''),width=8,height=4)
G
dev.off()

#!/usr/bin/Rscript
args=commandArgs(trailingOnly = TRUE)
library(data.table)
library(tidyverse)
library(gridExtra)

setwd("/project2/nobrega/grace/metabolic_traits_GWAS/")
rsid<-args[1]
trait=args[2]
cmd=paste("awk -F ' ' '$2==\"",rsid,"\" {print $1}' /project2/nobrega/grace/LD/GTEx_v8/*.bim",sep='')
chrom=as.numeric(system(cmd,intern=TRUE))

#Create data
cmd=paste("plink --bfile /project2/nobrega/grace/LD/GTEx_v8_females/",chrom," --ld-snp ",rsid," --ld-window 200000 --ld-window-r2 0 --maf 0.01 --out F_LD --r2",sep='')
system(cmd)  
cmd=paste("plink --bfile ~/midway/LD/GTEx_v8_males/",chrom," --ld-snp ",rsid," --ld-window 200000 --ld-window-r2 0 --maf 0.01 --out M_LD --r2",sep='')
system(cmd)  

#Load data
#Load LD
ld_F<-fread("F_LD.ld")
snp_pos<-ld_F$BP_A[1]
ld_F$rsid<-ld_F$SNP_B
ld_F$pos<-ld_F$BP_B
ld_F$chr<-ld_F$CHR_B
ld_F<-ld_F[,c("chr","pos","rsid","R2")]

ld_M<-fread("M_LD.ld")
ld_M$rsid<-ld_M$SNP_B
ld_M$pos<-ld_M$BP_B
ld_M$chr<-ld_M$CHR_B
ld_M<-ld_M[,c("chr","pos","rsid","R2")]

ld<-merge(ld_F,ld_M,by=c("chr","pos","rsid"),suffixes=c("_F","_M"))

#Load sumstats and find rsid name from location in hg19 
pos<-fread(paste("/project2/nobrega/grace/LD/GTEx_v7_LD/v7.GTEx.eur.",chrom,".bim",sep=''))
colnames(pos)<-c("chr","rsid","zero","pos","alt","ref")
pos<-pos[,c("rsid","pos")]
pos$pos<-as.character(pos$pos)

sumstats_F=fread(cmd=paste("zcat /scratch/midway2/gthansen/metabolic_traits_GWAS/",trait,".female.tsv.bgz",sep=''))
sumstats_F$chr<-sapply(strsplit(sumstats_F$variant,':'),'[[',1)
sumstats_F=sumstats_F[sumstats_F$chr==chrom,]
sumstats_F$pos<-sapply(strsplit(sumstats_F$variant,':'),'[[',2)
sumstats_F<-merge(pos,sumstats_F,by="pos")
gc()

sumstats_M=fread(cmd=paste("zcat /scratch/midway2/gthansen/metabolic_traits_GWAS/",trait,".male.tsv.bgz",sep=''))
sumstats_M$chr<-sapply(strsplit(sumstats_M$variant,':'),'[[',1)
sumstats_M=sumstats_M[sumstats_M$chr==chrom,]
sumstats_M$pos<-sapply(strsplit(sumstats_M$variant,':'),'[[',2)
sumstats_M<-merge(pos,sumstats_M,by="pos")
sumstats<-merge(sumstats_F,sumstats_M,by=c("pos","rsid"),suffixes=c("_F","_M"))
sumstats$pos<-NULL
gc()

dat<-merge(ld,sumstats,by="rsid")
#dat<-merge(dat,WHR_sex,by="rsid",suffixes=c("",paste("_",sex,sep='')))
dat$LD_F=character()
dat$LD_F[dat$R2_F<=0.2]="0-0.2"
dat$LD_F[dat$R2_F>0.2 & dat$R2_F<=0.4]="0.2-0.4"
dat$LD_F[dat$R2_F>0.4 & dat$R2_F<=0.6]="0.4-0.6"
dat$LD_F[dat$R2_F>0.6 & dat$R2_F<=0.8]="0.6-0.8"
dat$LD_F[dat$R2_F>0.8]="0.8-1.0"

dat$LD_M=character()
dat$LD_M[dat$R2_M<=0.2]="0-0.2"
dat$LD_M[dat$R2_M>0.2 & dat$R2_M<=0.4]="0.2-0.4"
dat$LD_M[dat$R2_M>0.4 & dat$R2_M<=0.6]="0.4-0.6"
dat$LD_M[dat$R2_M>0.6 & dat$R2_M<=0.8]="0.6-0.8"
dat$LD_M[dat$R2_M>0.8]="0.8-1.0"

dat<-dat[abs(dat$pos-snp_pos)<=200000,]
dat$pos=dat$pos/1000000

#Plot data
F<-ggplot(dat)+
  geom_point(aes(x=pos,y=-log10(pval_F),color=LD_F),size=2)+
  scale_y_continuous(limits=c(0,18))+
  scale_x_continuous(name="Position (MB)")+
  scale_color_manual(limits=c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"),values = c("dodgerblue4","dodgerblue1","green2","goldenrod2","firebrick3"))+
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme_minimal()+
  theme(panel.border = element_blank(), axis.line.x = element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22))

M<-ggplot(dat)+
  geom_point(aes(x=pos,y=-log10(pval_M),color=LD_M),size=2)+
  scale_y_continuous(limits=c(0,18))+
  scale_x_continuous(name="Position (MB)")+
  scale_color_manual(limits=c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"),values = c("dodgerblue4","dodgerblue1","green2","goldenrod2","firebrick3"))+
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme_minimal()+
  theme(panel.border = element_blank(), axis.line.x = element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22))

write.table(dat,paste("locus_zooms/",trait,"_",rsid,"_sex_locuszoom.txt",sep=''),quote=FALSE,sep='\t',row.names=FALSE)

pdf(paste("~/midway/metabolic_traits_GWAS/locus_zooms/",trait,"_",rsid,"_sex_locuszoom.pdf",sep=''),width=10,height=10)
grid.arrange(F,M,ncol=1)
dev.off()


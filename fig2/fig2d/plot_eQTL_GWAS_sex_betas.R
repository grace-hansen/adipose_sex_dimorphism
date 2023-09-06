#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(eulerr)
library(yarrr) #For color funs
library(optparse)
library(gridExtra)
trait<-"WHR"
tissue<-"adipose_subcutaneous"
data_source<-"GTEx_v8"
sex<-"F"
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")
setwd("~/medusa/papers/TWAS/TWAS/")

######### Author: Grace Hansen #########
#This script plots the beta values of variants in sex-specific eQTL analyses and sex-specific GWAS.

############################ For color manipulation ############################
darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
################################################################################

#Grab significant genes
genes<-fread(paste("~/midway/expression/",trait,"_",sex,"/",tissue,"/results/posthoc/sig_genes",sep=''))

#Grab top eQTL variant for each gene
TWAS<-fread(paste("~/midway/expression/",trait,"_",sex,"/",tissue,"/results/GTEx_v8.all.dat.top",sep=''))

#Grab p-value values of eQTL for each gene in each sex
get_vals<-function() {
  eQTL_B_F<-numeric(nrow(TWAS))
  eQTL_B_M<-numeric(nrow(TWAS))
  for (i in 1:nrow(TWAS)) {
    chr<-TWAS$CHR[i]
    gene<-TWAS$ID[i]
    var<-TWAS$EQTL.ID[i]
    B_F<-system(paste("awk -F ' ' '$1==\"",gene,"\"&&$2==\"",var,"\" {print $0}' ~/midway/QTL_analyses/eQTL/GTEx_F/results/chr",chr,"_eQTL.txt | cut -f4 -d' '",sep=''),intern=TRUE)
    B_M<-system(paste("awk -F ' ' '$1==\"",gene,"\"&&$2==\"",var,"\" {print $0}' ~/midway/QTL_analyses/eQTL/GTEx_M/results/chr",chr,"_eQTL.txt | cut -f4 -d' '",sep=''),intern=TRUE)
    if (length(B_F)>0 && length(B_M)>0) {
      eQTL_B_F[i]<-B_F
      eQTL_B_M[i]<-B_M
    }
  }
  #Grab p-values of eQTL variant for GWAS in each sex
  GWAS_B_F<-numeric(nrow(TWAS))
  GWAS_B_M<-numeric(nrow(TWAS))
  for (i in 1:nrow(TWAS)) {
    chr<-TWAS$CHR[i]
    gene<-TWAS$ID[i]
    var<-TWAS$EQTL.ID[i]
    B_F<-system(paste("zcat ~/midway/expression/WHR_F/Whradjbmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only.txt.gz | grep ",var,": | cut -f9 -d' '",sep=''),intern=TRUE)
    B_M<-system(paste("zcat ~/midway/expression/WHR_M/Whradjbmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only.txt.gz | grep ",var,": | cut -f9 -d' '",sep=''),intern=TRUE)
    if (length(B_F)>0 && length(B_M)>0) {
      GWAS_B_F[i]<-B_F
      GWAS_B_M[i]<-B_M
    }
  }
  return(list(eQTL_B_F,eQTL_B_M,GWAS_B_F,GWAS_B_M))
}
#vals<-get_vals()
#eQTL_B_F<-vals[[1]]
#eQTL_B_M<-vals[[2]]
#GWAS_B_F<-vals[[3]]
#GWAS_B_M<-vals[[4]]


#Load betas so I don't have to find them again
plot_dat<-fread("eQTL_GWAS_betas.txt")

#Normalize GWAS betas by h2
F_h2<-as.numeric(system("grep Observed ~/midway/ldsc/WHR_F/WHR_F_h2.log | cut -f5 -d' '",intern=TRUE))
M_h2<-as.numeric(system("grep Observed ~/midway/ldsc/WHR_M/WHR_M_h2.log | cut -f5 -d' '",intern=TRUE))
plot_dat$GWAS_B_F<-plot_dat$GWAS_B_F/F_h2
plot_dat$GWAS_B_M<-plot_dat$GWAS_B_M/M_h2
  
# Plot data as scatterplots: one for eQTL, one for GWAS, with regression line and line with slope=1

#plot_dat<-cbind(TWAS[,c("ID","EQTL.ID")],eQTL_B_F,eQTL_B_M,GWAS_B_F,GWAS_B_M)
#write.table(plot_dat,"eQTL_GWAS_betas.txt",quote=FALSE,row.names=FALSE,sep='\t')
plot_dat<-as.data.frame(plot_dat[rowSums(plot_dat!=0)==6,])
for (i in 3:6) {
  plot_dat[,i]<-as.numeric(plot_dat[,i])
}

eF<-abs(plot_dat$eQTL_B_F)
eM<-abs(plot_dat$eQTL_B_M)
slope_E<-lm(eM ~ eF)$coefficients[2]
int_E<-lm(eM ~ eF)$coefficients[1]
E<-ggplot(plot_dat,aes(x=abs(eQTL_B_F),y=abs(eQTL_B_M)))+
  geom_abline(linetype="dotted")+
  geom_abline(slope=slope_E,intercept=int_E,size=1)+
  geom_point(size=2,color=darken(rgb(pony_colors[3,1:3])))+
  annotate("text",x=0.3,y=1,label=paste("slope=",round(slope_E,digits=4),sep=''),size=5)+
  scale_x_continuous("abs(eQTL beta), females",limits=c(0,1.2))+
  scale_y_continuous("abs(eQTL beta), males",limits=c(0,1.2))+
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))

gF<-abs(plot_dat$GWAS_B_F)
gM<-abs(plot_dat$GWAS_B_M)
slope_G<-lm(gM ~ gF)$coefficients[2]
int_G<-lm(gM ~ gF)$coefficients[1]
G<-ggplot(plot_dat,aes(x=abs(GWAS_B_F),y=abs(GWAS_B_M)))+
  geom_abline(linetype="dotted")+
  geom_abline(slope=slope_G,intercept=int_G,size=1)+
  geom_point(size=2,color=rgb(pony_colors[7,1:3]))+
  annotate("text",x=0.05,y=0.25,label=paste("slope=",round(slope_G,digits=4),sep=''),size=5)+
  scale_x_continuous("normalized abs(GWAS beta), females",limits=c(0,0.35))+
  scale_y_continuous("normalized abs(GWAS beta), males",limits=c(0,0.35))+
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16))

pdf("eQTL_GWAS_sex_betas_scatter.pdf",width=4,height=8)
grid.arrange(E,G,nrow=2)
dev.off()

# Plot data as violin plots
eQTL<-pivot_longer(plot_dat,c(eQTL_B_F,eQTL_B_M),names_to="sex",values_to="eQTL")
GWAS<-pivot_longer(plot_dat,c(GWAS_B_F,GWAS_B_M),names_to="sex",values_to="GWAS")
eQTL$sex[eQTL$sex=="eQTL_B_F"]<-"F"
eQTL$sex[eQTL$sex=="eQTL_B_M"]<-"M"
GWAS$sex[GWAS$sex=="GWAS_B_F"]<-"F"
GWAS$sex[GWAS$sex=="GWAS_B_M"]<-"M"
plot_dat<-eQTL
plot_dat$GWAS<-GWAS$GWAS
plot_dat$GWAS_B_F<-NULL
plot_dat$GWAS_B_M<-NULL

E<-ggplot(plot_dat,aes(x=sex,y=abs(eQTL),group=EQTL.ID))+
  geom_violin(aes(x=sex,y=abs(eQTL),group=sex),fill=rgb(pony_colors[3,1:3]),draw_quantiles=c(0.5))+
  geom_point()+
  geom_line(alpha=0.2)+
  scale_y_continuous("abs(eQTL beta)",limits=c(0,1.5))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))
E

G<-ggplot(plot_dat,aes(x=sex,y=abs(GWAS),group=EQTL.ID))+
  geom_violin(aes(x=sex,y=abs(GWAS),group=sex),fill=rgb(pony_colors[7,1:3]),draw_quantiles=c(0.5))+
  geom_point()+
  geom_line(alpha=0.2)+
  scale_y_continuous("abs(GWAS beta)",limits=c(0,0.06))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))

grid.arrange(E,G,nrow=1)

pdf("eQTL_GWAS_sex_betas_violin.pdf",width=8,height=6)
grid.arrange(E,G,nrow=1)
dev.off()

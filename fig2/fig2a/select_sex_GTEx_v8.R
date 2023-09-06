library(data.table)
all<-fread("~/midway/expression/WHR/adipose_subcutaneous/adipose_subcutaneous_subject_IDs",header=FALSE)

#All phenotypes, GTEx v8
phenos<-fread("~/midway/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

#Subjects with adipose_subcutaneous gene exp and genotype data
phenos<-phenos[phenos$SUBJID %in% all$V1]

#Is this all the subject
if (length(unique(phenos$SUBJID))==nrow(all)) {
  males=phenos[phenos$SEX==1,'SUBJID']$SUBJID
  females=phenos[phenos$SEX==2,'SUBJID']$SUBJID
  males_to_exclude<-sample(males,length(males)-length(females),replace=FALSE)
  write(males,"~/midway/GTEx_v8_males",sep='\n')
  write(females,"~/midway/GTEx_v8_females",sep='\n')
  write(c(females,males_to_exclude),"~/midway/GTEx_v8_females_extra_males",sep='\n')
} else {
  print("Some samples are missing")
}

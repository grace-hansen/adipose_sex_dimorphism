#1/usr/bin/R
  
BPs<-c("A","T","G","C")

# Generate random barcodes
barcodes=character()
for (i in 1:200000) {
  barcode<-paste(sample(BPs,10,replace=TRUE),collapse='')
  barcodes[i]<-barcode
}

#Make sure barcodes are unique
barcodes<-unique(barcodes)

#Make sure barcodes don't contain Kpn1 or Xba1 binding sites
barcodes<-barcodes[!(grepl("TCTAGA",barcodes))]
barcodes<-barcodes[!(grepl("GGTACC",barcodes))]

#Make sure barcodes don't contain 4+ same base pair
barcodes<-barcodes[!(grepl("AAAA",barcodes))]
barcodes<-barcodes[!(grepl("TTTT",barcodes))]
barcodes<-barcodes[!(grepl("GGGG",barcodes))]
barcodes<-barcodes[!(grepl("CCCC",barcodes))]

#Make sure barcodes aren't an 'xxxyyy' structure
for (i in BPs) {
  for (j in BPs) {
    seq<-paste(c(i,i,i,j,j,j),collapse='')
    barcodes<-barcodes[!(grepl(seq,barcodes))]
  }
}

#Make sure barcodes don't end in TCT, TCTA, TCTAG
TCT<-logical(length=length(barcodes))
for (i in 1:length(barcodes)) {
  if (substr(barcodes[i],8,10)!="TCT") {
    TCT[i]<-TRUE
  }
  if (substr(barcodes[i],7,10)=="TCTA") {
    TCT[i]<-FALSE
  }
  if (substr(barcodes[i],6,10)=="TCTAG") {
    TCT[i]<-FALSE
  }
}
barcodes<-barcodes[TCT]

#Make sure barcodes don't end in GGT,GGTA,GGTAC
GGT<-logical(length=length(barcodes))
for (i in 1:length(barcodes)) {
  if (substr(barcodes[i],8,10)!="GGT") {
    GGT[i]<-TRUE
  }
  if (substr(barcodes[i],7,10)=="GGTA") {
    GGT[i]<-FALSE
  }
  if (substr(barcodes[i],6,10)=="GGTAC") {
    GGT[i]<-FALSE
  }
}
barcodes<-barcodes[GGT]

write(barcodes,"~/projects/MPRA/pseudorandom_barcodes.txt",sep='\n')

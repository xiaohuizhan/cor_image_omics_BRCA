dat<-read.table("Compare_Large_Area.txt",header=T,sep="\t")
library(ggplot2)
library(plyr)

dat$Name=factor(dat$Name, levels=as.character(dat$Name))
ggplot(dat, aes(x=Name,y=-log(dat$FDR.B.H,base=10),fill=Type))+
  geom_bar(stat="identity")+
  ylab("-log(FDR)")+ 
  scale_fill_manual(values=c("#e71010", "#3378ff"))+
  theme(axis.text.x = element_text(angle =90, size=6,hjust=0.5, vjust = 0.5))
  

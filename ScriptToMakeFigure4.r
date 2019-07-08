

setwd("~/Documents/Image/image_Jun/CPTAC_final/data_analysis/revise_manuscript/re-analysis/Enrichment_analysis_based_on_protein/positive_correlation/")

library(ggplot2)
library(tidyr)
dat<-read.table("enriched_BP_summary.txt",sep="\t",header=T)
dat$image=factor(dat$image,levels=c("Small_Area","Large_Area","Small_Major_Axis","Large_Major_Axis","Small_Minor_Axis","Large_Minor_Axis","Small_Aspect_Ratio","Large_Aspect_Ratio","Small_Mean_Distance","Large_Mean_Distance","Small_Max_Distance","Large_Max_Distance","Small_Min_Distance","Large_Min_Distance"))

name<-read.table("BP_name_rank.txt",header=F,sep="\t")
dat$Name=factor(dat$Name, levels=as.character(name$V1))

pdf("enrichment.pdf")

ggplot(dat, aes(x=image, y=Name,color=FDR.B.H)) + 
#geom_point(aes(size=100*FDR.B.H)) + 
geom_point() + 
scale_size_area() +scale_colour_gradient(low="#AB82FF",high="#EE2C2C")+
labs(title='The Biological Process Enrichment')+
theme(axis.title.x =element_text(size = 10),axis.title.y = element_text(size =10),plot.title = element_text(hjust =0.5,size=14))+ 
theme(axis.text.x = element_text(angle = 90, size=6,hjust=0.5, vjust = 0.5))+
theme(axis.text.y = element_text(size=6))

dev.off()
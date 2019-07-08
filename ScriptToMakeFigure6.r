
setwd("~/Documents/Image/image_Jun/CPTAC_final/data_analysis/revise_manuscript/re-analysis/survival/2019.7.7/")

source("https://bioconductor.org/biocLite.R")
biocLite("forestplot")
biocLite("haven")
library("forestplot")
library("haven")

Forestplot<-read.table("prognosis_image_feature.txt",header=T,sep="\t")

lower=round(Forestplot$cox_lower, 2) 
upper=round(Forestplot$cox_upper, 2) 
CI<-paste(lower,upper,sep="-")
CI<-paste("(",CI,")",sep="")
hr=round(Forestplot$cox_hr, 2)
HR95CI<-paste(hr,CI,sep="")

dat<-data.frame(Forestplot$Morphological_feature,HR95CI,hr,lower,upper)
colnames(dat)<-c("Varianles","HR95CI","HR","LowerCI","UpperCI")
attach(dat)

tabletext<-cbind(
  c("Varianles",as.character(dat$Varianles)),
  c("HR95CI", as.character(dat$HR95CI)))

forestplot(as.matrix(dat[,1:3]),HR,LowerCI,UpperCI,graph.pos=2,zero=1,
	graphwidth=unit(50,"mm"),lineheight="auto",boxsize=0.1,
	xticks=(c(0,0.5,1.0,1.5,2.0,2.5)),col=fpColors(all.elements="black"),
	title="Forestplot of image feature")




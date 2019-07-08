library(survival)
library(ggplot2)
library(survminer)
rm(list=ls())
setwd("~/Documents/Image/image_Jun/CPTAC_final/data_analysis/revise_manuscript/re-analysis/survival/")

clinical<-read.table("BRCA_clinicalMatrix",header=T,sep="\t")
newclinical<-clinical[which(clinical$OS.time != "NA"), ] 
new<-data.frame(newclinical$sampleID,newclinical$OS.time,newclinical$OS)
cli<-c()
for(i in 1:dim(new)[1]){
     tmp<-strsplit(as.character(new$newclinical.sampleID[i]),"-")
     if(tmp[[1]][4]=="01"){
         samp<-paste(tmp[[1]][1],tmp[[1]][2],tmp[[1]][3],sep=".")	
         cli<-rbind(cli,c(samp,as.character(as.matrix(new[i,-1]))))
     }
}
cli<-as.data.frame(cli)
colnames(cli)<-c("patient","OS.time","OS")

image<-read.table("image_BRCA.txt",sep="\t",header=T)
new.image<-as.matrix(image)
i<-1
res<-c()
while (i< 100) {
     tmp<-new.image[c(i,i+1),]
     dat<-apply(tmp, 2, sum)
     res<-rbind(res,dat)
     i<-i+2
 }
fin<-rbind(res,new.image[101:150,])
name<-read.table("imFeasName.txt",header=F,sep="\t")
rownames(fin)<-as.character(name$V1)
image<-fin

new.image<-t(image)
patient<-rownames(new.image)
new.image<-data.frame(patient,new.image)

res<-merge(cli,new.image,by="patient")

low<-floor(0.2*dim(res)[1])  #20th percentiles 向下取整
high<-ceiling(0.8*dim(res)[1])#80th percentiles 向下取整


##image.analysis<-c("Area5","DistMin1","Minor5","DistMean1", "DistMax5")


selected.image<-c("Area1","Area5","Major1","Major5","Minor1","Minor5","Ratio1","Ratio5","DistMean1","DistMean5",
"DistMax1", "DistMax5","DistMin1","DistMin5")

rename.image<-c("Small_Nucleus_Area","Large_Nucleus_Area","Small_Major_Axis","Large_Major_Axis","Small_Minor_Axis","Large_Minor_Axis",
"Small_Aspect_Ratio","Large_Aspect_Ratio","Small_Mean_Distance","Large_Mean_Distance",
"Small_Max_Distance", "Large_Max_Distance","Small_Min_Distance","Large_Min_Distance")

out<-c()
for(i in 1:length(selected.image)){
     
     value<-sort(as.numeric(as.character(res[,which(colnames(res)==selected.image[i])])))
     select<-c()
    
     for(j in low:high ){
         cutoff<-value[j]
         classification  <- ifelse(res[,which(colnames(res)==selected.image[i])] > cutoff,1,0)
         dat<-data.frame(as.numeric(as.character(res$OS.time)),as.numeric(as.character(res$OS)),classification)
         colnames(dat)<-c("time","status","classification")
         sdf<-survdiff(Surv(dat$time,dat$status) ~ classification)
         p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
         select<-rbind(select,c(cutoff,j/dim(res)[1],p.val))
     }
    select<-as.data.frame(select)
    colnames(select)<-c("cutoff","percent","pvalue")
    select<-na.omit(select)

    min.pvalue<-min(as.numeric(as.character(as.matrix(select$pvalue))))
    num<-which(as.numeric(as.character(as.matrix(select$pvalue)))==min.pvalue)
    fin.cutoff<-select$cutoff[num]
   
    ## cox regression

    classify <- ifelse(res[,which(colnames(res)==selected.image[i])] > fin.cutoff, 1,0)
    fin.dat<-data.frame(as.numeric(as.character(res$OS.time)),as.numeric(as.character(res$OS)),classify)
    colnames(fin.dat)<-c("time","status","classify")
    
    cox<-coxph(Surv(fin.dat$time,fin.dat$status) ~ classify)
    hr <- as.numeric(summary(cox)$coefficients[2])
    lower <- as.numeric(summary(cox)$conf.int[3])
    upper <- as.numeric(summary(cox)$conf.int[4])

    ## output
    out<-rbind(out,c(rename.image[i],as.matrix(select[num,]), hr,lower,upper)) 

    
    ## line plot

    pva<--log(min.pvalue,10)
    point<-data.frame(fin.cutoff,pva)

    pdf(paste("~/Documents/Image/image_Jun/CPTAC_final/data_analysis/revise_manuscript/re-analysis/survival/line_chart_2/",rename.image[i],".pdf",sep=""))
    
    print(ggplot(data = select, aes(x = cutoff, y = -log(pvalue,10))) +
    geom_line(color="#92c5de")+
    geom_point(data=point,aes(x=fin.cutoff,y=pva,color="#c51b7d"),size=2,shape=20)+
    xlab("Image feature value")+
    ylab("-log10pvalue")+ ggtitle(rename.image[i])+
    theme(plot.title = element_text(hjust = 0.5))+guides(color=FALSE))
    dev.off()
    
     ## survival 

     filename<-paste("~/Documents/Image/image_Jun/CPTAC_final/data_analysis/revise_manuscript/re-analysis/survival/survival_2/",rename.image[i],".pdf",sep="")
     

    fin.dat$classify[fin.dat$classify==1]<-"high_image_value"
    fin.dat$classify[fin.dat$classify==0]<-"low_image_value"
     gg<-ggsurvplot(survfit(Surv(fin.dat$time, fin.dat$status)~classify, data=fin.dat),
     pval=min.pvalue,pval.size=3,pval.coord=c(0.2,0.05),legend=c(0.8,0.8),xlab="Time(days)",
     title=rename.image[i],ggtheme = theme_minimal() + 
     theme(plot.title = element_text(hjust = 0.5, face = "bold")))
     ggsave(filename,print(gg))


}

out<-as.data.frame(out)
colnames(out)<-c("Morphological_feature","cutoff","percent","log_rank_pvalue","cox_hr","cox_lower","cox_upper")
write.table(out,"survival_analysis_output.txt",quote=F,sep="\t",row.names=F)


## feature selection: log_rank_test pvalue<0.001 ; HR >1.2 or <1/1.2
tmp<-out[which(as.numeric(as.character(out$log_rank_pvalue))<0.001),]

unfavourite<-tmp[which(as.numeric(as.character(tmp$cox_hr))>1.2),]
favourite<-tmp[which(as.numeric(as.character(tmp$cox_hr))<1/1.2),]

prog<-rbind(unfavourite,favourite)
write.table(prog,"prognosis_image_feature.txt",quote=F,sep="\t",row.names=F)


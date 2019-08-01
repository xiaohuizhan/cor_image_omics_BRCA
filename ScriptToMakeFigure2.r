
setwd("~/Documents/Image/image_Jun/CPTAC_final/data_analysis/revise_manuscript/re-analysis/RPPA/")
RPPA<-read.table("BRCA.rppa.txt",header=T,sep="\t")

## deal with sampleID && select tumor sample
name<-colnames(RPPA)
patient_ID<-c("Composite.Element.REF")
site<-c(1)
for(i in 2:length(name)){
    tmp<-strsplit(as.character(name[i]),"\\.")
    tmp1<-strsplit(as.character(tmp[[1]][4]),"")
    tmp2<-paste(tmp1[[1]][1],tmp1[[1]][2],sep="")
        if(tmp2<10){
         samp<-paste(tmp[[1]][1],tmp[[1]][2],tmp[[1]][3],sep=".")	
         patient_ID<-c(patient_ID,samp)
         site<-c(site,i)
    }
}
new.rppa<-RPPA[,site]
colnames(new.rppa)<-patient_ID


## select patient with matched image data
image<-read.table("image.txt",header=T,sep="\t")
patient<-colnames(image)
rppa<-new.rppa[,patient]
rppa<-data.frame(new.rppa[,1],rppa)

## proteins with missing values in more than 20% of the samples were excluded 

missdata<-c()
for(i in 1:dim(rppa)[1]){
     tmp<-sum(is.na(rppa[i,]))
     if(tmp>(ceiling(dim(rppa)[2]*0.2)+1)){
      missdata<-c(missdata,i)
     }
}
rppa<-rppa[-missdata,]

##proteins with missing values in less than 20% of samples were imputed 

library("mice")
tempData <- mice(rppa,m=5,meth='pmm') 
summary(tempData) 
res.rppa <- complete(tempData) 
anyNA(res.rppa) 
rownames(res.rppa)<-rownames(rppa)
colnames(res.rppa)[1]<-"ID"
write.table(res.rppa,"rppa.txt",quote=F,sep="\t")


## correlation calculation for RPPA data

cor<-c("rppa_ID","image","rppa.cor")
cor<-data.frame(cor)
cor<-t(cor)
write.table(cor,"rppa_cor.txt",sep="\t",quote=F,row.names=F,col.names=F)

for(i in 1:dim(res.rppa)[1]){
    x<-as.numeric(as.character(as.matrix(res.rppa[i,-1])))
    for(j in 1:dim(image)[1]){
         z<-as.numeric(as.character(image[j,]))
         cor.res.rppa<-cor.test(x, z, alternative = "two.sided", method = "spearman",conf.level = 0.95)
         out<-t(c(as.character(res.rppa[i,1]),rownames(image)[j],
              as.numeric(as.character(cor.res.rppa$estimate))))
        write.table(out,"rppa_cor.txt",sep="\t",quote=F,row.names=F,col.names=F,append=T)
   }   
}

##rppa
rppa<-read.table("rppa_cor.txt",sep="\t",header=T)
gene<-sapply(strsplit(as.character(rppa$rppa_ID),'\\|'), "[", 1)
Antibody<-sapply(strsplit(as.character(rppa$rppa_ID),'\\|'), "[", 2)
new.rppa<-data.frame(gene,Antibody,rppa)

## cptac
cptac<-read.csv("correlation.txt",sep="\t",header=T)

### compare the correlation of common gene between rppa and cptac
 
select.image<-c("Area1","Area5","Major1","Major5","Minor1","Minor5","Ratio1","Ratio5","DistMean1","DistMean5","DistMax1","DistMax5","DistMin1","DistMin5")
title.name<-c("Small_Area","Large_Area","Small_Major_Axis","Large_Major_Axis","Small_Minor_Axis","Large_Minor_Axis","Small_Aspect_Ratio","Large_Aspect_Ratio","Small_Mean_Distance","Large_Mean_Distance","Small_Max_Distance","Large_Max_Distance","Small_Min_Distance","Large_Min_Distance")

library(ggplot2)

for(i in 1:length(select.image)){

    x<-new.rppa[which(as.character(new.rppa$image)==select.image[i]),]
    y<-cptac[which(as.character(cptac$image)==select.image[i]),]

    # deal with the gene has multi-antibody for RPPA data
    m<-data.frame(table(x$gene)) 
    n<-m[order(m[,2],decreasing=T),]
    k<-n[which(as.numeric(as.character(n$Freq))>1),]
    w<-n[which(as.numeric(as.character(n$Freq))==1),]
    colnames(w)[1]<-"gene"
    tmp.rppa<-merge(w,x,by="gene")[,-2]

    for(j in 1:dim(k)[1]){
    
        tmp<-x[which(as.character(x$gene)==as.character(k$Var1[j])),]
        cor<-as.numeric(as.character(tmp$rppa.cor))
        if(min(cor)>0){
        tmp.rppa<-rbind(tmp.rppa,tmp[which(cor==max(cor)),])
        }
        if(max(cor)<0){
        tmp.rppa<-rbind(tmp.rppa,tmp[which(cor==min(cor)),])
        }
    }

    ## compare the correlation between rppa-image and cptac-image
    z<-as.data.frame(intersect(tmp.rppa$gene,y$gene))
    colnames(z)<-"gene"
    fil.rppa<-merge(z,tmp.rppa,by="gene")
    fil.cptac<-merge(z,y,by="gene")
    dat<-merge(fil.rppa[,c(1,5)],fil.cptac[,c(1,5)],by="gene")

     output<-ggplot(dat, aes(x=rppa.cor, y=protein.cor)) + 
                   geom_point(shape = 16, size = 3, color="red",show.legend = FALSE, alpha = .5)+
                   stat_smooth(method="lm")+ 
                   scale_y_continuous(limits=c(-0.5, 0.5))+
                   scale_x_continuous(limits=c(-0.5, 0.5))+
                   labs(x = "", y = "")+
                   ggtitle(title.name[i])+
                   xlab("ρ_RPPA")+
                   ylab("ρ_CPTAC")+ 
                   theme(plot.title = element_text(hjust = 0.5))+
                   guides(color=FALSE)
     ggsave(output,filename=paste(title.name[i],".pdf",sep=""))

}




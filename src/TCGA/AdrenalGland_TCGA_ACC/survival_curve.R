library(survsim)
library(ggplot2)
library(rlang)

#-------------read parameter file----------------                                        
fp<-read.table("data/Rparameter.txt",sep="\t")
d<-as.matrix(fp)
devideType=d[2,2]
highCut=NULL
lowCut=NULL
if(devideType=="user"){
    lowCut=as.numeric(d[3,2])
    highCut=as.numeric(d[4,2])
}
subtype=d[5,2]

data<-read.table("data/FPKM_ALL_sort_geneName_clinic_subtype_head_data.txt",sep="\t")
dm<-as.matrix(data)
m<-dm[,5:ncol(dm)]
m<-m[,which(m[2,]=="Primary Tumor")]
gene_type<-NULL
if(subtype!="All"){
    m<-m[,which(m[3,]==subtype)]
}

if(dm[nrow(dm),1]=="GeneName"){
    png("result/TCGA_high_low_expression_survival.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()

    png("result/TCGA_high_low_expression.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()
}

rowMax <- nrow(dm)
geneName=dm[rowMax,1]

if(devideType=="user"){
    mL<-m[,which(as.numeric(m[rowMax,])<=lowCut)]
    mH<-m[,which(as.numeric(m[rowMax,])>highCut)]
    m<-cbind(mL,mH)
    gene_type<-ifelse(as.numeric(m[rowMax,])<=lowCut,0,1)
}else if(devideType=="median"){
    gene_median<-median(as.numeric(m[rowMax,]))
    gene_type<-ifelse(as.numeric(m[rowMax,])<=gene_median,0,1)
    mL<-m[,which(as.numeric(m[rowMax,])<=gene_median)]
    mH<-m[,which(as.numeric(m[rowMax,])>gene_median)]
}else if(devideType=="tertile"){
    temp<-sort(as.numeric(m[rowMax,]))
    lowCut<-temp[(length(temp)/3.0)]
    highCut<-temp[(length(temp)*2/3.0)]
    print(lowCut)
    print(highCut)
    mL<-m[,which(as.numeric(m[rowMax,])<=lowCut)]
    mH<-m[,which(as.numeric(m[rowMax,])>highCut)]
    m<-cbind(mL,mH)
    gene_type<-ifelse(as.numeric(m[rowMax,])<=lowCut,0,1)
}else if(devideType=="quartile"){
    temp<-sort(as.numeric(m[rowMax,]))
    lowCut<-temp[(length(temp)/4.0)]
    highCut<-temp[(length(temp)*3/4.0)]
    mL<-m[,which(as.numeric(m[rowMax,])<=lowCut)]
    mH<-m[,which(as.numeric(m[rowMax,])>highCut)]
    m<-cbind(mL,mH)
    gene_type<-ifelse(as.numeric(m[rowMax,])<=lowCut,0,1)
}else if(devideType=="mean_sd"){
    m_8 <- as.numeric(m[rowMax,])
    meanVal <- mean(log2(m_8+0.01))
    sd <- sd(log2(m_8+0.01))
    mL<-m[,which(log2(m_8+0.01)<=(meanVal-sd))]
    mH<-m[,which(log2(m_8+0.01)>(meanVal+sd))]
    m<-cbind(mL,mH)
    gene_type<-ifelse(as.numeric(log2(as.numeric(m[rowMax,])+0.01))<=(meanVal-sd),0,1)    
}


                                        #------------------------draw curve----------------------------------
if(ncol(mH)>0 && ncol(mL)>0){
    vital<-m[3,]
    vital[vital=="Alive"]<-0
    vital[vital=="Dead"]<-1
    event<-ifelse(vital%in%c("0","FALSE"),FALSE,TRUE)
    followupMonth<-(as.numeric(m[4,])/30.0)

    sfit=Surv(as.numeric(followupMonth),as.numeric(vital))
    fit=survfit(sfit~gene_type)

    sdf<-survdiff(sfit~factor(gene_type))
    p.val<-1-pchisq(sdf$chisq,length(sdf$n)-1)
    print("type:")
    print(subtype)
    print("p-value:")
    print(p.val)
                                        #-----------------------------------------draw survival
    png("result/TCGA_high_low_expression_survival.png",width=6,height=5,res=100,units='in')
    par(mar=c(5,6,1,2),mgp=c(3.5,1,0))
    plot(fit,col=c(1:2),frame=F,mark.time=T,lwd=3,cex.axis=2.5,yaxt="n",xlab="Months",ylab="Survival Probability",cex.lab=2.5,mgp=c(3.5,1.5,0))
    legend(0,0.3,legend=c(paste("High ",ncol(mH)),paste("Low  ",ncol(mL))),col=c(2:1),bty="n",lwd=4,cex=1.5)
    legend(0,0.12, legend=paste(expression("logRank p="),format(p.val,digits=2)) ,bty="n",cex=1.5)
    axis(1,cex.axis=2,lwd=3,labels=NA)
    axis(2,cex.axis=2.5,at=seq(0,1,by=0.5),labels=c("0","0.5","1.0"),lwd=3)
    dev.off()
    
                                        #-----------------------------------------draw expression
    c1 <- data.frame(group="High",value=as.numeric(mH[rowMax,]))
    c2 <- data.frame(group="Low",value=as.numeric(mL[rowMax,]))
    df <- rbind(c1,c2)
    ggplot(df,aes(x=group, y=value,fill=group,colour=group))+stat_boxplot(geom="errorbar",width=0.4,lwd=1,color="black")+geom_boxplot(width=0.6,lwd=0.5,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.2),size=2,alpha=1)+labs(y=paste(geneName," (FPKM)"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=25,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_blank(),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(plot.margin=unit(c(0.2,0.2,0.2,0.5),"in"))+geom_hline(yintercept=1,linetype="dashed",size=1,color="red")+scale_colour_manual(values=c("red","black"))+scale_fill_manual(values=c("gray","gray"))
    ggsave("result/TCGA_high_low_expression.png",width=6,height=5,units='in')

                                        #-------------------------write table----------------------------
    mH<-t(mH)
    mH<-mH[order(mH[,3]),]
    mH<-cbind(1:nrow(mH),mH[,1],mH[,2],mH[,3],mH[,4],mH[,rowMax])
    mH <- rbind(c("index","TCGA","Type","Vital","Followup","FPKM"),mH)
    write.table(mH,"result/TCGA_high_low_expression_H.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    mL<-t(mL)
    mL<-mL[order(mL[,3]),]
    mL<-cbind(1:nrow(mL),mL[,1],mL[,2],mL[,3],mL[,4],mL[,rowMax])
    mL <- rbind(c("index","TCGA","Type","Vital","Followup","FPKM"),mL)
    write.table(mL,"result/TCGA_high_low_expression_L.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

}else{
    png("result/TCGA_high_low_expression_survival.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()

    png("result/TCGA_high_low_expression.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()
    write.table(NULL,"result/TCGA_high_low_expression_H.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(NULL,"result/TCGA_high_low_expression_L.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

                                        #-------------------------subtypes distribution figure----------------------------
for(i in (6:rowMax-1)){
    typeName <- dm[i,4]
    typeName <- gsub("\\/","_",typeName)
    typeName <- gsub("\\ ","_",typeName)
    typeName <- gsub("\\=","_",typeName)
    tmp <- dm[,5:ncol(dm)]
    tmp <- tmp[c(i,1,rowMax),which(tmp[2,]!="Solid Tissue Normal")]
    rownames(tmp) <- c("group","barcode","value")
    colnames(tmp) <- NULL
    df <- data.frame(t(tmp))
    write.table(df,paste("result/TCGA_boxplot_tumor_",typeName,".data",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    df$value <- log2(as.numeric(tmp[3,])+0.01)
    factorSize = max(nchar(as.character(df$group)),na.rm=TRUE)
    factorSize1 = nchar(as.character(df$group))[1]
    print(factorSize1)
    angle=0
    hjust=0.5
    leftM=0.4
    if(factorSize>4 && factorSize<12){
        angle=30
        hjust=1
    }else if(factorSize>=12){
        angle=60
        hjust=1
        leftM=1
    }

    ggplot(df,aes(x=group, y=value,fill=group))+stat_boxplot(geom="errorbar",width=0.4)+geom_boxplot(lwd=1.0,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.2),size=2,alpha=1)+labs(y=paste(geneName," log2(FPKM)"),x=typeName)+theme(panel.grid.major = element_blank(), panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=25,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_text(face="plain",size=25,vjust=0.4),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(axis.text.x=element_text(angle=angle,hjust=hjust))+theme(plot.margin=unit(c(0.4,0.1,0.1,leftM),"in"))
    ggsave(paste("result/TCGA_boxplot_tumor_",typeName,".png",sep=""),height=5+factorSize/6,width=6,units='in')


                                        #-------------------------wilcoxon test----------------------------
    types <- unique(tmp[1,])
    types <- types[ !types %in% NA ]
    wil <- NULL
    firstRow <- c("pvalue")
    for(j in (1:length(types))){
        firstRow <- c(firstRow,types[j])
        nextRow <- c(types[j])
        for(k in (1:length(types))){
            d1 <- as.numeric(tmp[3,which(tmp[1,]==types[j])])
            d2 <- as.numeric(tmp[3,which(tmp[1,]==types[k])])
            val<-round(as.numeric(wilcox.test(d1,d2,paired=FALSE)$p.value),digits=3)
            nextRow <- c(nextRow,val)
        }
        wil <- rbind(wil,nextRow)
    }
    wil <- rbind(firstRow,wil)
    write.table(wil,paste("result/TCGA_boxplot_tumor_",typeName,".wilcox",sep=""),quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE)
}

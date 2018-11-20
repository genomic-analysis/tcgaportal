library(survsim)
library(ggplot2)
library(rlang)

#-------------read parameter file----------------                                        
fp<-read.table("Rparameter.txt",sep="\t")
d<-as.matrix(fp)
devideType=d[2,2]
highCut=NULL
lowCut=NULL
if(devideType=="user"){
    lowCut=as.numeric(d[3,2])
    highCut=as.numeric(d[4,2])
}
subtype=d[5,2]

data<-read.table("data/TCGA_clinical_expression.txt",sep="\t")
dm<-as.matrix(data)
m<-dm[,5:ncol(dm)]
m<-m[,which(m[2,]=="Primary Tumor")]
MGAM2_type<-NULL
if(subtype!="All"){
    m<-m[,which(m[3,]==subtype)]
}

if(nrow(dm)<7){
    png("TCGA_survival_curve.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()

    png("TCGA_high_low_expression.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()
}

geneName=dm[7,1]

if(devideType=="user"){
    mL<-m[,which(as.numeric(m[7,])<=lowCut)]
    mH<-m[,which(as.numeric(m[7,])>highCut)]
    m<-cbind(mL,mH)
    MGAM2_type<-ifelse(as.numeric(m[7,])<=lowCut,0,1)
}else if(devideType=="median"){
    MGAM2_median<-median(as.numeric(m[7,]))
    MGAM2_type<-ifelse(as.numeric(m[7,])<=MGAM2_median,0,1)
    mL<-m[,which(as.numeric(m[7,])<=MGAM2_median)]
    mH<-m[,which(as.numeric(m[7,])>MGAM2_median)]
}else if(devideType=="tertile"){
    temp<-sort(as.numeric(m[7,]))
    lowCut<-temp[(length(temp)/3.0)]
    highCut<-temp[(length(temp)*2/3.0)]
    print(lowCut)
    print(highCut)
    mL<-m[,which(as.numeric(m[7,])<=lowCut)]
    mH<-m[,which(as.numeric(m[7,])>highCut)]
    m<-cbind(mL,mH)
    MGAM2_type<-ifelse(as.numeric(m[7,])<=lowCut,0,1)
}else if(devideType=="quartile"){
    temp<-sort(as.numeric(m[7,]))
    lowCut<-temp[(length(temp)/4.0)]
    highCut<-temp[(length(temp)*3/4.0)]
    mL<-m[,which(as.numeric(m[7,])<=lowCut)]
    mH<-m[,which(as.numeric(m[7,])>highCut)]
    m<-cbind(mL,mH)
    MGAM2_type<-ifelse(as.numeric(m[7,])<=lowCut,0,1)
}else if(devideType=="mean_sd"){
    m_8 <- as.numeric(m[7,])
    meanVal <- mean(log2(m_8+0.01))
    sd <- sd(log2(m_8+0.01))
    mL<-m[,which(log2(m_8+0.01)<=(meanVal-sd))]
    mH<-m[,which(log2(m_8+0.01)>(meanVal+sd))]
    m<-cbind(mL,mH)
    MGAM2_type<-ifelse(as.numeric(log2(as.numeric(m[7,])+0.01))<=(meanVal-sd),0,1)    
}

                                        #------------------------draw curve----------------------------------
if(ncol(mH)>0 && ncol(mL)>0){
    vital<-m[5,]
    vital[vital=="Alive"]<-0
    vital[vital=="Dead"]<-1
    event<-ifelse(vital%in%c("0","FALSE"),FALSE,TRUE)
    followupMonth<-(as.numeric(m[6,])/30.0)

    sfit=Surv(as.numeric(followupMonth),as.numeric(vital))
    fit=survfit(sfit~MGAM2_type)

    sdf<-survdiff(sfit~factor(MGAM2_type))
    p.val<-1-pchisq(sdf$chisq,length(sdf$n)-1)
    print("type:")
    print(subtype)
    print("p-value:")
    print(p.val)
    png("TCGA_survival_curve.png",width=6,height=5,res=100,units='in')
    par(mar=c(5,6,1,2),mgp=c(3.5,1,0))
    plot(fit,col=c(1:2),frame=F,mark.time=T,lwd=3,cex.axis=2.5,yaxt="n",xlab="Months",ylab="Survival Probability",cex.lab=2.5,mgp=c(3.5,1.5,0))
    legend(0,0.4,legend=c("high","low"),col=c(2:1),bty="n",lwd=4,cex=2)
    legend(0,0.17, legend=paste(expression("logRank p="),format(p.val,digits=2)) ,bty="n",cex=2)
    axis(1,cex.axis=2,lwd=3,labels=NA)
    axis(2,cex.axis=2.5,at=seq(0,1,by=0.5),labels=c("0","0.5","1.0"),lwd=3)
    dev.off()
    
                                        #-----------------------------------------draw expression
    c1 <- data.frame(group="High",value=as.numeric(mH[7,]))
    c2 <- data.frame(group="Low",value=as.numeric(mL[7,]))
    df <- rbind(c1,c2)
    ggplot(df,aes(x=group, y=value,fill=group,colour=group))+stat_boxplot(geom="errorbar",width=0.4,lwd=1,color="black")+geom_boxplot(width=0.6,lwd=0.5,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.2),size=2,alpha=1)+labs(y=paste(geneName," (FPKM)"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=25,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_blank(),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(plot.margin=unit(c(0.2,0.2,0.2,0.5),"in"))+geom_hline(yintercept=1,linetype="dashed",size=1,color="red")+scale_colour_manual(values=c("red","black"))+scale_fill_manual(values=c("gray","gray"))
    ggsave("TCGA_high_low_expression.png",width=6,height=5,units='in')

                                        #-------------------------write table----------------------------
    mHT<-t(mH)
    mHT<-mHT[order(mHT[,3]),]
    mHT<-cbind(1:nrow(mHT),mHT[,1],mHT[,2],mHT[,4],mHT[,6],mHT[,7])

    write.table(mHT,"data/TCGA_expression_H.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    mLT<-t(mL)
    mLT<-mLT[order(mLT[,3]),]
    mLT<-cbind(1:nrow(mLT),mLT[,1],mLT[,2],mLT[,4],mLT[,6],mLT[,7])

    write.table(mLT,"data/TCGA_expression_L.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    
}else{
    png("TCGA_survival_curve.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()

    png("TCGA_high_low_expression.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()

    png("TCGA_primary_tumor_subtype_boxplot.png",width=6,height=5,res=100,units='in')
    plot(0)
    dev.off()

}

                                        #-------------------------draw subtype distribution boxplot and save wilcoxo-test pvalue----------------------------
tmp <- dm[,5:ncol(dm)]
noNormalM <- tmp[c(3,1,7),which(tmp[3,]!="Solid Tissue Normal")]
rownames(noNormalM) <- c("group","barcode","value")
colnames(noNormalM) <- NULL
df <- data.frame(t(noNormalM))
print(geneName)
df$value <- log2(as.numeric(noNormalM[3,]))
ggplot(df,aes(x=group, y=value,fill=group))+stat_boxplot(geom="errorbar",width=0.4)+geom_boxplot(lwd=1.0,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.2),size=2,alpha=1)+labs(y=paste(geneName," log2(FPKM)"))+theme(panel.grid.major = element_blank(), panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=25,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_blank(),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(axis.text.x=element_text(angle=30,hjust=1))+theme(plot.margin=unit(c(0.1,0.1,0.1,0.4),"in"))
ggsave("TCGA_primary_tumor_subtype_boxplot.png",height=6,width=6,units='in')

m1 <- as.numeric(noNormalM[3,which(noNormalM[1,]=="prox.-inflam")])
m2 <- as.numeric(noNormalM[3,which(noNormalM[1,]=="prox.-prolif.")])
m3 <- as.numeric(noNormalM[3,which(noNormalM[1,]=="TRU")])
a1<-round(as.numeric(wilcox.test(m1,m1,paired=FALSE)$p.value),digits=3)
a2<-round(as.numeric(wilcox.test(m1,m2,paired=FALSE)$p.value),digits=3)
a3<-round(as.numeric(wilcox.test(m1,m3,paired=FALSE)$p.value),digits=3)
b1<-round(as.numeric(wilcox.test(m2,m1,paired=FALSE)$p.value),digits=3)
b2<-round(as.numeric(wilcox.test(m2,m2,paired=FALSE)$p.value),digits=3)
b3<-round(as.numeric(wilcox.test(m2,m3,paired=FALSE)$p.value),digits=3)
c1<-round(as.numeric(wilcox.test(m3,m1,paired=FALSE)$p.value),digits=3)
c2<-round(as.numeric(wilcox.test(m3,m2,paired=FALSE)$p.value),digits=3)
c3<-round(as.numeric(wilcox.test(m3,m3,paired=FALSE)$p.value),digits=3)
d <- rbind(c("pval","prox.-inflam","prox.-prolif.","TRU"),
c("prox.-inflam",a1,a2,a3),
c("prox.-prolif.",b1,b2,b3),
c("TRU",c1,c2,c3))
write.table(d,"TCGA_primary_tumor_subtype_boxplot.wilcox",quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE)


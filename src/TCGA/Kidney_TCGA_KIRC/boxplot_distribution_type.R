library(ggplot2)
library(rlang)

f=read.table("TCGA_clinical_expression.txt",sep="\t")
dm<-as.matrix(f)

#############
temp <- t(dm)
write.table(temp,"temp",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
#############

if(nrow(dm)<9){
    png("boxplot_distribution_stage.png",height=5,width=5,res=100,units='in')
    plot(0)
    dev.off()
}

geneName <- dm[9,1]
#dm<-dm[,which(dm[2,]=="Primary Tumor")]

m1F=cbind(dm[9,which(dm[2,]=="Primary Tumor")])
m2F=cbind(dm[9,which(dm[2,]=="Solid Tissue Normal")])

a1<-data.frame(group="Primary Tumor",value=log2(as.numeric(m1F)))
b1<-data.frame(group="Solid Tissue Normal",value=log2(as.numeric(m2F)))
df <-rbind(a1,b1)

png("boxplot_distribution_type.png",height=6,width=5,res=100,units='in')
ggplot(df,aes(x=group, y=value,fill=group))+stat_boxplot(geom="errorbar",width=0.5)+geom_boxplot(lwd=0.5,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.2),size=2,color="black",alpha=1)+labs(x="Sample Type",y=paste(geneName," (log2(FPKM))"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=25,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_text(face="plain",size=25,vjust=0.4),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(axis.text.x=element_text(angle=30,hjust=0.95))
dev.off()


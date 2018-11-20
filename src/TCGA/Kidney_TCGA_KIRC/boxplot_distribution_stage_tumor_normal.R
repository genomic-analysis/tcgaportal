library(ggplot2)
library(rlang)

f=read.table("TCGA_clinical_expression.txt",sep="\t")
dm<-as.matrix(f)

if(nrow(dm)<9){
    png("boxplot_distribution_stage.png",height=5,width=5,res=100,units='in')
    plot(0)
    dev.off()
}

geneName <- dm[9,1]
tm<-dm[,which(dm[2,]=="Primary Tumor")]
nm<-dm[,which(dm[2,]=="Solid Tissue Normal")]

t1=cbind(tm[9,which(tm[7,]=="Stage I")])
t2=cbind(tm[9,which(tm[7,]=="Stage II")])
t3=cbind(tm[9,which(tm[7,]=="Stage III")])
t4=cbind(tm[9,which(tm[7,]=="Stage IV")])
a1<-data.frame(group="Stage I",Type="tumor",value=log2(as.numeric(t1)))
b1<-data.frame(group="Stage II",Type="tumor",value=log2(as.numeric(t2)))
c1<-data.frame(group="Stage III",Type="tumor",value=log2(as.numeric(t3)))
d1<-data.frame(group="Stage IV",Type="tumor",value=log2(as.numeric(t4)))

n1=cbind(nm[9,which(nm[7,]=="Stage I")])
n2=cbind(nm[9,which(nm[7,]=="Stage II")])
n3=cbind(nm[9,which(nm[7,]=="Stage III")])
n4=cbind(nm[9,which(nm[7,]=="Stage IV")])
a2<-data.frame(group="Stage I",Type="normal",value=log2(as.numeric(n1)))
b2<-data.frame(group="Stage II",Type="normal",value=log2(as.numeric(n2)))
c2<-data.frame(group="Stage III",Type="normal",value=log2(as.numeric(n3)))
d2<-data.frame(group="Stage IV",Type="normal",value=log2(as.numeric(n4)))

df <-rbind(a1,a2,b1,b2,c1,c2,d1,d2)

png("boxplot_distribution_stage_tumor_normal.png",height=6,width=6,res=100,units='in')
ggplot(df,aes(x=group, y=value,fill=Type))+stat_boxplot(geom="errorbar", lwd=1.5,width=0.45,position=position_dodge(width=0.7))+geom_boxplot(width=0.5,lwd=1.0,position=position_dodge(width=0.7),outlier.shape=NA)+geom_jitter(position=position_dodge(width=0.7),size=2,color="black",alpha=1)+labs(y=paste(geneName," log2(FPKM)"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1.5,colour="black"),axis.text.x=element_text(size=28,face="plain",color="black"),axis.text.y=element_text(size=28,face="plain",color="black"),axis.title.x=element_blank(),axis.title.y=element_text(face="plain",size=28,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(axis.text.x=element_text(angle=30,hjust=1))+theme(plot.margin=unit(c(0.2,0.2,0.2,0.4),"in"))+scale_fill_manual(values = c("#FF0000","#008000","#FF0000","#008000","#FF0000","#008000","#FF0000","#008000")) + theme(legend.position = "bottom",legend.text=element_text(size=15))
dev.off()

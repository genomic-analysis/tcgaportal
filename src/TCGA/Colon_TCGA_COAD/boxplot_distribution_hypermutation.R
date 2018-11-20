library(ggplot2)
library(rlang)

f=read.table("TCGA_clinical_expression.txt",sep="\t")
dm<-as.matrix(f)

if(nrow(dm)<9){
    png("boxplot_distribution_hypermutation.png",height=5,width=5,res=100,units='in')
    plot(0)
    dev.off()
}

geneName <- dm[6,1]
dm<-dm[,which(dm[2,]=="Primary Tumor")]

m0F=cbind(dm[6,which(dm[3,]=="0")])
m1F=cbind(dm[6,which(dm[3,]=="1")])

a1<-data.frame(group="Non-hypermutated",value=log2(as.numeric(m0F)))
b1<-data.frame(group="Hypermutated",value=log2(as.numeric(m1F)))
df <-rbind(a1,b1)

#Wilcoxon ranking test                                        
m0<-as.numeric(m0F)
m1<-as.numeric(m1F)

pval<-round(as.numeric(wilcox.test(m0,m1,paired=FALSE)$p.value),3)
png("boxplot_distribution_hypermutation.png",height=5,width=5,res=100,units='in')
ggplot(df,aes(x=group, y=value,fill=group))+stat_boxplot(geom="errorbar",width=0.5)+geom_boxplot(lwd=0.5,width=0.5,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.15),size=2,color="black",alpha=1)+labs(x="",y=paste(geneName," (log2(FPKM))"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=16,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_text(face="plain",size=25,vjust=0.4),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(axis.text.x=element_text(angle=0,hjust=0.5))+theme(plot.margin=unit(c(0.2,0.2,0.1,0.1),"in"))+annotate("text", x = -Inf, y = Inf, label = paste(" Wilcoxo test pval:",pval), hjust = 0, vjust = 1, parse = FALSE)
           
dev.off()

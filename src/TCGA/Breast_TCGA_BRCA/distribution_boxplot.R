library(ggplot2)
library(rlang)

f=read.table("TCGA_clinical_expression.txt",sep="\t")
dm<-as.matrix(f)
geneName <- dm[6,1]
dm<-dm[,which(dm[3,]=="Primary Tumor")]
LumaF=cbind(dm[6,which(dm[2,]=="LumA")])
LumbF=cbind(dm[6,which(dm[2,]=="LumB")])
BasalF=cbind(dm[6,which(dm[2,]=="Basal")])
Her2F=cbind(dm[6,which(dm[2,]=="Her2")])
NormalF=cbind(dm[6,which(dm[2,]=="Normal")])

a1<-data.frame(group="LumA",value=log2(as.numeric(LumaF)))
b1<-data.frame(group="LumB",value=log2(as.numeric(LumbF)))
c1<-data.frame(group="Basal",value=log2(as.numeric(BasalF)))
d1<-data.frame(group="Her2",value=log2(as.numeric(Her2F)))
e1<-data.frame(group="Normal-like",value=log2(as.numeric(NormalF)))
df <-rbind(a1,b1,c1,d1,e1)

LumA = length(as.numeric(LumaF))
LumA_H = sum(as.numeric(LumaF)>1.0)
LumB = length(as.numeric(LumbF))
LumB_H = sum(as.numeric(LumbF)>1.0)
Basal = length(as.numeric(BasalF))
Basal_H = sum(as.numeric(BasalF)>1.0)
Her2 = length(as.numeric(Her2F))
Her2_H = sum(as.numeric(Her2F)>1.0)
Normal = length(as.numeric(NormalF))
Normal_H = sum(as.numeric(NormalF)>1.0)
sum_H=LumA_H+LumB_H+Basal_H+Her2_H+Normal_H
sum = length(LumaF)+length(LumbF)+length(BasalF)+length(Her2F)+length(NormalF)
LumA_matrix=matrix(c(LumA_H,sum_H,LumA,sum),2,2)
LumB_matrix=matrix(c(LumB_H,sum_H,LumB,sum),2,2)
Basal_matrix=matrix(c(Basal_H,sum_H,Basal,sum),2,2)
Her2_matrix=matrix(c(Her2_H,sum_H,Her2,sum),2,2)
Normal_matrix=matrix(c(Normal_H,sum_H,Normal,sum),2,2)

fisher_LumA=paste("p=",format(fisher.test(LumA_matrix,alternative="two.sided")$p.value,digits=2),sep="")
fisher_LumB=paste("p=",format(fisher.test(LumB_matrix)$p.value,digits=2),sep="")
fisher_Basal=paste("p=",format(fisher.test(Basal_matrix)$p.value,digits=2),sep="")
fisher_Her2=paste("p=",format(fisher.test(Her2_matrix)$p.value,digits=2))
fisher_Normal=paste("p=",format(fisher.test(Normal_matrix)$p.value,digits=2),sep="")

print(fisher_LumA)
print(fisher_LumB)
print(fisher_Basal)
print(fisher_Her2)
print(fisher_Normal)


#Wilcoxon ranking test
LumA<-as.numeric(LumaF)
LumB<-as.numeric(LumbF)
Basal<-as.numeric(BasalF)
Her2<-as.numeric(Her2F)
Normal<-as.numeric(NormalF)
a1<-as.numeric(wilcox.test(LumA,c(LumB,Basal,Her2,Normal),paired=FALSE)$statistic)
a2<-as.numeric(wilcox.test(LumB,c(LumA,Basal,Her2,Normal),paired=FALSE)$statistic)
a3<-as.numeric(wilcox.test(Basal,c(LumA,LumB,Her2,Normal),paired=FALSE)$statistic)
a4<-as.numeric(wilcox.test(Her2,c(LumA,LumB,Basal,Normal),paired=FALSE)$statistic)
a5<-as.numeric(wilcox.test(Normal,c(LumA,LumB,Basal,Her2),paired=FALSE)$statistic)

##
a1<-round(as.numeric(wilcox.test(LumA,LumA,paired=FALSE)$p.value),digits=3)
a2<-round(as.numeric(wilcox.test(LumA,LumB,paired=FALSE)$p.value),digits=3)
a3<-round(as.numeric(wilcox.test(LumA,Basal,paired=FALSE)$p.value),digits=3)
a4<-round(as.numeric(wilcox.test(LumA,Her2,paired=FALSE)$p.value),digits=3)
a5<-round(as.numeric(wilcox.test(LumA,Normal,paired=FALSE)$p.value),digits=3)

b1<-round(as.numeric(wilcox.test(LumB,LumA,paired=FALSE)$p.value),digits=3)
b2<-round(as.numeric(wilcox.test(LumB,LumB,paired=FALSE)$p.value),digits=3)
b3<-round(as.numeric(wilcox.test(LumB,Basal,paired=FALSE)$p.value),digits=3)
b4<-round(as.numeric(wilcox.test(LumB,Her2,paired=FALSE)$p.value),digits=3)
b5<-round(as.numeric(wilcox.test(LumB,Normal,paired=FALSE)$p.value),digits=3)

c1<-round(as.numeric(wilcox.test(Basal,LumA,paired=FALSE)$p.value),digits=3)
c2<-round(as.numeric(wilcox.test(Basal,LumB,paired=FALSE)$p.value),digits=3)
c3<-round(as.numeric(wilcox.test(Basal,Basal,paired=FALSE)$p.value),digits=3)
c4<-round(as.numeric(wilcox.test(Basal,Her2,paired=FALSE)$p.value),digits=3)
c5<-round(as.numeric(wilcox.test(Basal,Normal,paired=FALSE)$p.value),digits=3)

d1<-round(as.numeric(wilcox.test(Her2,LumA,paired=FALSE)$p.value),digits=3)
d2<-round(as.numeric(wilcox.test(Her2,LumB,paired=FALSE)$p.value),digits=3)
d3<-round(as.numeric(wilcox.test(Her2,Basal,paired=FALSE)$p.value),digits=3)
d4<-round(as.numeric(wilcox.test(Her2,Her2,paired=FALSE)$p.value),digits=3)
d5<-round(as.numeric(wilcox.test(Her2,Normal,paired=FALSE)$p.value),digits=3)

e1<-round(as.numeric(wilcox.test(Normal,LumA,paired=FALSE)$p.value),digits=3)
e2<-round(as.numeric(wilcox.test(Normal,LumB,paired=FALSE)$p.value),digits=3)
e3<-round(as.numeric(wilcox.test(Normal,Basal,paired=FALSE)$p.value),digits=3)
e4<-round(as.numeric(wilcox.test(Normal,Her2,paired=FALSE)$p.value),digits=3)
e5<-round(as.numeric(wilcox.test(Normal,Normal,paired=FALSE)$p.value),digits=3)

#print(c(a1,a2,a3,a4,a5))


print(c("",   "LumA","LumB","Basal","Her2","Normal"),quote=FALSE)
print(c("LumA",   a1,   a2,     a3,    a4,     a5),quote=FALSE)
print(c("LumB",   b1,   b2,     b3,    b4,     b4),quote=FALSE)
print(c("Basal",  c1,   c2,     c3,    c4,     c5),quote=FALSE)
print(c("Her2",   d1,   d2,     d3,    d4,     d5),quote=FALSE)
print(c("Normal", e1,   e2,     e3,    e4,     e5),quote=FALSE)

wil <- rbind(c("pval",  "LumA","LumB","Basal","Her2","Normal"),c("LumA",   a1,   a2,     a3,    a4,     a5),c("LumB",   b1,   b2,     b3,    b4,     b4),c("Basal",  c1,   c2,     c3,    c4,     c5),c("Her2",   d1,   d2,     d3,    d4,     d5),c("Normal", e1,   e2,     e3,    e4,     e5))
write.table(wil,paste("TCGA_tumor_sample_boxplot.wilcox",sep=""),quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE)


png("TCGA_tumor_sample_boxplot.png",height=6,width=6,res=100,units='in')
ggplot(df,aes(x=group, y=value,fill=group))+stat_boxplot(geom="errorbar",width=0.5)+geom_boxplot(lwd=0.5,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.2),size=2,color="black",alpha=1)+labs(x="MGAM2",y=paste(geneName," log2(FPKM)"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=25,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_blank(),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(axis.text.x=element_text(angle=30,hjust=1))
dev.off()

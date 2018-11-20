library(ggplot2)
library(rlang)

f=read.table("TCGA_clinical_expression.txt",sep="\t")
dm<-as.matrix(f)

if(nrow(dm)<9){
    png("boxplot_distribution_mRNA_cluster.png",height=5,width=5,res=100,units='in')
    plot(0)
    dev.off()
}

geneName <- dm[9,1]
dm<-dm[,which(dm[2,]=="Primary Tumor")]

m1F=cbind(dm[9,which(dm[3,]=="1")])
m2F=cbind(dm[9,which(dm[3,]=="2")])
m3F=cbind(dm[9,which(dm[3,]=="3")])
m4F=cbind(dm[9,which(dm[3,]=="4")])

a1<-data.frame(group="m1",value=log2(as.numeric(m1F)))
b1<-data.frame(group="m2",value=log2(as.numeric(m2F)))
c1<-data.frame(group="m3",value=log2(as.numeric(m3F)))
d1<-data.frame(group="m4",value=log2(as.numeric(m4F)))
df <-rbind(a1,b1,c1,d1)

m1 = length(as.numeric(m1F))
m1_H = sum(as.numeric(m1F)>1.0)
m2 = length(as.numeric(m2F))
m2_H = sum(as.numeric(m2F)>1.0)
m3 = length(as.numeric(m3F))
m3_H = sum(as.numeric(m3F)>1.0)
m4 = length(as.numeric(m4F))
m4_H = sum(as.numeric(m4F)>1.0)

sum_H=m1_H+m2_H+m3_H+m4_H
sum = length(m1F)+length(m2F)+length(m3F)+length(m4F)
m1_matrix=matrix(c(m1_H,sum_H,m1,sum),2,2)
m2_matrix=matrix(c(m2_H,sum_H,m2,sum),2,2)
m3_matrix=matrix(c(m3_H,sum_H,m3,sum),2,2)
m4_matrix=matrix(c(m4_H,sum_H,m4,sum),2,2)

fisher_m1=paste("p=",format(fisher.test(m1_matrix,alternative="two.sided")$p.value,digits=2),sep="")
fisher_m2=paste("p=",format(fisher.test(m2_matrix)$p.value,digits=2),sep="")
fisher_m3=paste("p=",format(fisher.test(m3_matrix)$p.value,digits=2),sep="")
fisher_m4=paste("p=",format(fisher.test(m4_matrix)$p.value,digits=2))

print(fisher_m1)
print(fisher_m2)
print(fisher_m3)
print(fisher_m4)

#Wilcoxon ranking test
m1<-as.numeric(m1F)
m2<-as.numeric(m2F)
m3<-as.numeric(m3F)
m4<-as.numeric(m4F)

a1<-as.numeric(wilcox.test(m1,c(m2,m3,m4),paired=FALSE)$statistic)
a2<-as.numeric(wilcox.test(m2,c(m1,m3,m4),paired=FALSE)$statistic)
a3<-as.numeric(wilcox.test(m3,c(m1,m2,m4),paired=FALSE)$statistic)
a4<-as.numeric(wilcox.test(m4,c(m1,m2,m3),paired=FALSE)$statistic)


##
a1<-round(as.numeric(wilcox.test(m1,m1,paired=FALSE)$p.value),digits=3)
a2<-round(as.numeric(wilcox.test(m1,m2,paired=FALSE)$p.value),digits=3)
a3<-round(as.numeric(wilcox.test(m1,m3,paired=FALSE)$p.value),digits=3)
a4<-round(as.numeric(wilcox.test(m1,m4,paired=FALSE)$p.value),digits=3)

b1<-round(as.numeric(wilcox.test(m2,m1,paired=FALSE)$p.value),digits=3)
b2<-round(as.numeric(wilcox.test(m2,m2,paired=FALSE)$p.value),digits=3)
b3<-round(as.numeric(wilcox.test(m2,m3,paired=FALSE)$p.value),digits=3)
b4<-round(as.numeric(wilcox.test(m2,m4,paired=FALSE)$p.value),digits=3)

c1<-round(as.numeric(wilcox.test(m3,m1,paired=FALSE)$p.value),digits=3)
c2<-round(as.numeric(wilcox.test(m3,m2,paired=FALSE)$p.value),digits=3)
c3<-round(as.numeric(wilcox.test(m3,m3,paired=FALSE)$p.value),digits=3)
c4<-round(as.numeric(wilcox.test(m3,m4,paired=FALSE)$p.value),digits=3)

d1<-round(as.numeric(wilcox.test(m4,m1,paired=FALSE)$p.value),digits=3)
d2<-round(as.numeric(wilcox.test(m4,m2,paired=FALSE)$p.value),digits=3)
d3<-round(as.numeric(wilcox.test(m4,m3,paired=FALSE)$p.value),digits=3)
d4<-round(as.numeric(wilcox.test(m4,m4,paired=FALSE)$p.value),digits=3)

print(c("",   "m1","m2","m3","m4"),quote=FALSE)
print(c("m1",   a1,   a2,     a3,    a4),quote=FALSE)
print(c("m2",   b1,   b2,     b3,    b4),quote=FALSE)
print(c("m3",  c1,   c2,     c3,    c4),quote=FALSE)
print(c("m4",   d1,   d2,     d3,    d4),quote=FALSE)
cat(a1)
cat("\n")
cat(a2)
cat("\n")
cat(a3)
cat("\n")
cat(a4)
cat("\n")
cat(b1)
cat("\n")
cat(b2)
cat("\n")
cat(b3)
cat("\n")
cat(b4)
cat("\n")
cat(c1)
cat("\n")
cat(c2)
cat("\n")
cat(c3)
cat("\n")
cat(c4)
cat("\n")
cat(d1)
cat("\n")
cat(d2)
cat("\n")
cat(d3)
cat("\n")
cat(d4)
cat("\n")

png("boxplot_distribution_mRNA_cluster.png",height=5,width=5,res=100,units='in')
ggplot(df,aes(x=group, y=value,fill=group))+stat_boxplot(geom="errorbar",width=0.5)+geom_boxplot(lwd=0.5,color="black",outlier.shape=NA)+geom_jitter(position=position_jitter(width=0.2),size=2,color="black",alpha=1)+labs(x="mRNA class",y=paste(geneName," (log2(FPKM))"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black",size=1),axis.ticks.x=element_line(size=1,colour="black"),axis.ticks.y=element_line(size=1,colour="black"),axis.text.x=element_text(size=25,face="plain",color="black"),axis.text.y=element_text(size=25,face="plain",color="black"),axis.title.x=element_text(face="plain",size=25,vjust=0.4),axis.title.y=element_text(face="plain",size=25,vjust=0.4),legend.position="none",line=element_line(size=2,colour="black"))+theme(axis.text.x=element_text(angle=0,hjust=0.5))+theme(plot.margin=unit(c(0,0,0.01,0.1),"in"))
dev.off()

library(survsim)
library(ggplot2)
library(rlang)

f1 <- read.table("TCGA_expressionA.txt",sep="\t")
m1 <- as.matrix(f1)
f2 <- read.table("TCGA_expressionB.txt",sep="\t")
m2 <- as.matrix(f2)

if( nrow(m1)>5 &&  nrow(m2)>5  ){
    d1 <- as.numeric(m1[6,which(m1[3,]=="Primary Tumor")])
    g1 <- m1[6,1]
    
    d2 <- as.numeric(m2[6,which(m2[3,]=="Primary Tumor")])
    g2 <- m2[6,1]
    
    pearson=cor.test(d1,d2, method="pearson")
    spearman=cor.test(d1,d2, method="spearman")
    png("TCGA_correlation.png",width=6,height=6,res=300,units='in')
    par(mar=c(7,7,1,1),mgp=c(4.0,1.8,0))
    plot(d1,d2,xlab=g1,ylab=g2,pch=19,cex=0.2,cex.lab=2.5,cex.axis=2.5,lwd=10)
    box(lwd=3)
    axis(side=1,lwd=3,lab=NA)
    axis(side=2,lwd=3,lab=NA)
    title(sub=paste(expression("Pearson = "),format(pearson$estimate,digits=3)), adj=0, line=4, font=2)
    title(sub=paste(expression("Spearman = "),format(spearman$estimate,digits=3)), adj=0, line=5, font=2)
    dev.off()
    print(pearson$estimate)
    print(spearman$estimate)
    
}else{
    png("TCGA_correlation.png")
    plot(0)
    dev.off()
}

library(survsim)

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

#------------- read survival data file----------------------
data<-read.table("METABRIC_clinical_expression.txt",sep="\t")
dm<-as.matrix(data)
m<-dm[,3:ncol(dm)]
MGAM2_type<-NULL
if(subtype!="All"){
m<-m[,which(m[2,]==subtype)]
}


if(devideType=="user"){
mL<-m[,which(m[7,]<lowCut)]
mH<-m[,which(m[7,]>highCut)]
m<-cbind(mL,mH)
MGAM2_type<-ifelse(as.numeric(m[7,])<lowCut,0,1)
}else if(devideType=="median"){
val<-m[7,]
val<-val[which(val!="NA")]
MGAM2_median<-median(as.numeric(val))
MGAM2_type<-ifelse(as.numeric(m[7,])<MGAM2_median,0,1)
mL<-m[,which(m[7,]<MGAM2_median)]
mH<-m[,which(m[7,]>MGAM2_median)]
}else if(devideType=="tertile"){
temp<-sort(m[7,])
lowCut<-temp[(length(temp)/3.0)]
highCut<-temp[(length(temp)*2/3.0)]
print(lowCut)
print(highCut)
mL<-m[,which(m[7,]<=lowCut)]
mH<-m[,which(m[7,]>highCut)]
m<-cbind(mL,mH)
MGAM2_type<-ifelse(as.numeric(m[7,])<lowCut,0,1)
}else if(devideType=="quartile"){
temp<-sort(m[7,])
lowCut<-temp[(length(temp)/4.0)]
highCut<-temp[(length(temp)*3/4.0)]
mL<-m[,which(m[7,]<=lowCut)]
mH<-m[,which(m[7,]>highCut)]
m<-cbind(mL,mH)
MGAM2_type<-ifelse(as.numeric(m[7,])<lowCut,0,1)
}


#------------------------draw curve----------------------------------
vital<-m[3,]
vital[vital=="LIVING"]<-0
vital[vital=="DECEASED"]<-1
event<-ifelse(vital%in%c("0","FALSE"),FALSE,TRUE)
followupMonth<-as.numeric(m[4,])

sfit=Surv(as.numeric(followupMonth),as.numeric(vital))
fit=survfit(sfit~MGAM2_type)

sdf<-survdiff(sfit~factor(MGAM2_type))
p.val<-1-pchisq(sdf$chisq,length(sdf$n)-1)
print("type:")
print(subtype)
print("p-value:")
print(p.val)
png("METABRIC_survival_curve.png",width=7,height=5,res=100,units='in')
par(mar=c(5,7,2,0),mgp=c(4,1,0))
plot(fit,col=c(1:2),xlim=c(0,500),frame=F,mark.time=T,lwd=3,yaxt="n",xaxt="n",xlab="months",ylab="Survival Probability",cex.lab=3)
legend(200,0.95,legend=c("high","low"),col=c(2:1),bty="n",lwd=4,cex=2)
legend(100, 1.1, legend=paste(expression("logRank p="),format(p.val,digits=2)) ,bty="n",cex=2)
axis(1,cex.axis=2,at=seq(0,400,by=200),labels=c("0","200","400"),cex.axis=3,lwd=3,mgp=c(0,2,0))
axis(2,cex.axis=2,at=seq(0,1,by=0.5),labels=c("0","0.5","1.0"),cex.axis=3,lwd=3)
dev.off()


#-------------------------write table----------------------------
mH<-t(mH)
mH<-cbind(1:nrow(mH),mH[,1],mH[,2],mH[,6],mH[,4],mH[,5])
mH<-mH[order(mH[,4]),]
write.table(mH,"METABRIC_expression_H.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
mL<-t(mL)
mL<-cbind(1:nrow(mL),mL[,1],mL[,2],mL[,6],mL[,4],mL[,5])
mL<-mL[order(mL[,4]),]
write.table(mL,"METABRIC_expression_L.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

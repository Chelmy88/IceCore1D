setwd("~/Documents/Paper/Ice_model/time_series")
library(splines)

EDCage<-as.matrix(read.table("EDC_AICC2012_official.csv",header=TRUE, sep="\t"))
EDCage<-EDCage[,c(1,2)]
EDCage[,1]<-3261-EDCage[,1]
EDCage[,2]<-EDCage[,2]+50
age_EDC<-as.data.frame(spline(EDCage,xout=seq(5,3250,by=5),method="natural"))
plot(EDCage)
lines(age_EDC,col="red")
write.table(age_EDC[,2], "EDC_age_forC.dat", sep = "\t", row.names = FALSE,col.names = FALSE)
EDCage[1:50,]


############
tempEDC<-as.matrix(read.table("LR04-EDC_temp_1Myr.dat",header=FALSE, sep="\t"))
accEDC<-as.matrix(read.table("LR04-EDC_acc_1Myr.dat",header=FALSE, sep="\t"))
par(mfrow=c(1,1))
temp_mod<-tempEDC
temp_mod[8900:9850]<-temp_mod[8900:9850]+2
temp_mod[8120:8600]<-temp_mod[8120:8600]+2

plot(temp_mod[8120:8600])
plot(tempEDC,type="l")
lines(temp_mod,col="blue")
plot(accEDC,type="l")

#####PRODUCE PLOTS####

tempCN<-as.matrix(read.table("export/temp_profile_test_EDC_CN.dat",header=FALSE, sep="\t"))
ageCN<-as.matrix(read.table("export/ageRel_test_EDC_CN.dat",header=FALSE, sep="\t"))
par(mfrow=c(1,2))
plot(age_EDC[,2],type="l")
lines(ageCN[,2],col="red")
plot(temp_EDC[,1]+273.15,type="l")
lines(rev(tempCN[1:3262]),col="red")



# 
# 
# 
EDC<-as.matrix(read.table("EDC temp profile.csv",header=FALSE, sep="\t"))
 EDC<-EDC[order(EDC[,2],decreasing=TRUE),]
 EDC0<-EDC
# 
tail(EDC0)
EDC0
EDC[1:61,]=EDC[1:61,]+0.08
EDCs1<-as.data.frame(spline(EDC[,c(2,1)],xout=seq(-7,-3194,by=-1),method="fmm"))
EDCs1<-EDCs1[,c(2,1)]
EDCs<-as.data.frame(spline(EDCs1[-c(c(50:85),94,c(1671:1695),c(828:850)),c(2,1)],xout=seq(-7,-3194,by=-1),method="fmm"))
EDCs<-EDCs[,c(2,1)]

par(mfrow=c(2,2))
plot(EDC,type="l",lwd=3,ylim=c(-3310,0))
plot(EDC[1:2000,],type="l",lwd=1)
lines(EDC0[1:80,],type="l",lwd=1,col="blue")
lines(EDCs[1:130,],col="red")

set.seed(1)
tail(EDCs)
y <- EDCs[1600:3188,1]
x <- EDCs[1600:3188,2]
## Fit a couple of quadratic splines with different degrees of freedom
f1 <- lm(y ~ bs(x, degree = 4,df=9))  # Defaults to 2 - 1 = 1 degree of freedom
x2 <- EDCs[3188,2]-c(0:80)
## Plot the splines
plot(EDCs[3000:3188,],type="l",xlim=c(-9,0),ylim=c(-3300,-3000))
lines(predict(f1, data.frame(x = x2)),x2, col = "red")
abline(v=-2.9,lty=2)
abline(h=EDCs[3188,2]-67,lty=2)
x3<-EDCs[3188,2]-c(1:67)
EDCf<-rbind(as.matrix(EDCs),as.matrix(cbind(predict(f1, data.frame(x = x3)),x3)))
plot(EDCf,type="l")#EDCf is from 7 to 3261


# 
# 
# 
# 
# tempCN
# 
# ageCN
# age_EDC
# 
# tempCN<-as.matrix(read.table("export/setest_temp_profile_CN.dat",header=FALSE, sep="\t"))
# tempCN<-cbind(rev(tempCN[1:3262]-273.15),c(0:-3261))
# plot(tempCN,type="l")
# lines(EDCf,col="blue")
# 
# ageCN<-as.matrix(read.table("export/setest_ageRel_CN.dat",header=FALSE, sep="\t"))
# 
# 
# EDCage[,2]<-EDCage[,2]+55
# lines(EDCage,col="red")
# EDCage<-as.data.frame(spline(EDCage,xout=seq(0,3192,by=1),method="natural"))
# 
# write.table((EDCage[,2]), "EDC_age_forC.dat", sep = "\t", row.names = FALSE,col.names = FALSE)
# write.table((rev(EDCf[,1]+273.15)), "EDC_temp_forC.dat", sep = "\t", row.names = FALSE,col.names = FALSE)
# 
# par(mfrow=c(1,2))
# plot(EDCf,type="l")
# temp_load_all<-as.matrix(read.table(sprintf("export/temp_profile_with_se_and_lat_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
# lines(rev(temp_load_all[1:3262])-273.15,c(-1:-3262),col=2)    
# temp_load_se<-as.matrix(read.table(sprintf("export/temp_profile_with_se_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
# lines(rev(temp_load_se[1:3262])-273.15,c(-1:-3262),col=3)  
# temp_load_lat<-as.matrix(read.table(sprintf("export/temp_profile_with_lat_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
# lines(rev(temp_load_lat[1:3262])-273.15,c(-1:-3262),col=4)  
# temp_load_without<-as.matrix(read.table(sprintf("export/temp_profile_without_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
# lines(rev(temp_load_without[1:3262])-273.15,c(-1:-3262),col=5)  
# temp_load_bump<-as.matrix(read.table(sprintf("export/temp_profile_test_bump_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
# lines(rev(temp_load_bump[1:3262])-273.15,c(-1:-3262),col=6)  
# 
# plot(-EDCf[,1]+(rev(temp_load_se[1:3255])-273.15))
# plot(-EDCf[,1]+(rev(temp_load_bump[1:3255])-273.15))
# plot(-EDCf[,1]+(rev(temp_load_lat[1:3255])-273.15))
# 
# plot(-temp_load_without+temp_load_se)
# plot(-temp_load_without+temp_load_bump)
# plot(-temp_load_without+temp_load_lat)
# lines(temp_load_without-temp_load_all)
# 
# 
# par(mar=c(4,4,3,1))
# par(mfrow=c(1,2))
# plot(EDCf,type="l",main="Profile",ylab="Depth (m)",xlab="Temperature (C)")
# temp_load_all<-as.matrix(read.table(sprintf("export/temp_profile_with_se_and_lat_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
# lines(rev(temp_load_all[1:3262])-273.15,c(-1:-3262),col=2)    
# legend("bottomleft",bty="n",legend=c("Measurment","Model"),col=c(1,2),lty=1,ncol=1)
# plot((EDCf[,2]),-EDCf[,1]+(rev(temp_load_se[1:3255])-273.15),type="l",main="Difference",xlab="Depth (m)",ylab="Temperature (C)",xlim=c(0,-3200))
# 
# 
# 
# 
# par(mfrow=c(1,1))
# th<-as.matrix(read.table(sprintf("LR04-EDC_thickness_1Myr.dat"),header=FALSE, sep="\t"))
# plot(1+abs(th[1:10000]-th[2:10001])/(31556926.*100.)/10^(-13),type="l")
# plot(1+abs(th[1:10000]-th[2:10001])/(31556926.*100.)/10^(-10),type="l")
# 
# 
# summary1<-as.matrix(read.table("export/Summary_CN1.dat",header=FALSE, sep="\t"))
# summary2<-as.matrix(read.table("export/Summary_CN.dat",header=FALSE, sep="\t"))
# 
# summary2<-summary2[1:45,1:4]
# 
# summary2
# colo = heat.colors(20)
# 
# par(mfrow=c(2,2))
# par(mar=c(2,2,1,1))
# plot(summary1[,1],log10(summary1[,3]),col=colo[100*summary1[,2]])
# text(summary1[,1], log10(summary1[,3]), labels=summary1[,2], cex= 0.7,pos=4)
# 
# plot(summary1[,2],log10(summary1[,3]),col=colo[10*summary1[,1]])
# text(summary1[,2], log10(summary1[,3]), labels=summary1[,1], cex= 0.7,pos=4)
# 
# plot(summary1[,1],(summary1[,4]),col=colo[100*summary1[,2]])
# text(summary1[,1], (summary1[,4]), labels=summary1[,2], cex= 0.7,pos=4)
# 
# plot(summary1[,2],(summary1[,4]),col=colo[10*summary1[,1]])
# text(summary1[,2], (summary1[,4]), labels=summary1[,1], cex= 0.7,pos=4)
# 
# par(mfrow=c(2,2))
# par(mar=c(2,2,1,1))
# 
# plot(summary2[,1],log10(summary2[,3]),col=colo[100*summary2[,2]])
# text(summary2[,1], log10(summary2[,3]), labels=summary2[,2], cex= 0.7,pos=4)
# 
# plot(summary2[,2],log10(summary2[,3]),col=colo[10*summary2[,1]])
# text(summary2[,2], log10(summary2[,3]), labels=summary2[,1], cex= 0.7,pos=4)
# 
# plot(summary2[,1],(summary2[,4]),col=colo[100*summary2[,2]])
# text(summary2[,1], (summary2[,4]), labels=summary2[,2], cex= 0.7,pos=4)
# 
# plot(summary2[,2],(summary2[,4]),col=colo[10*summary2[,1]])
# text(summary2[,2], (summary2[,4]), labels=summary2[,1], cex= 0.7,pos=4)
# 
# a=0.035
# print(a,digits=4)
# par(mfrow=c(1,1))
# 
# 
# Qlist<-seq(0.050,0.06,by=0.001)
# folder="Run2_k_P_v_fi"
# m=0.4
# 
# plot(EDCf,type="l",xlim=c(-15,0),ylim=c(-3300,-2500))
# for(QG in Qlist){
#   temp_load<-as.matrix(read.table(sprintf("export/%s/temp_profile_m_%.2f_QG_%2.1f_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
#   lines(rev(temp_load[1:3262])-273.15,c(-1:-3262),col=QG*900)    
# }
# legend("topright",legend=Qlist,col=Qlist*900,lty=1,ncol=2)
# lines(EDCf,type="l",lwd=2,col="red")
# 
# plot(EDCf,type="l",xlim=c(-50,-30),ylim=c(-1800,-700),lwd=2)
# for(QG in Qlist){
#   temp_load<-as.matrix(read.table(sprintf("export/%s/temp_profile_m_%.2f_QG_%2.1f_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
#   lines(rev(temp_load[1:3262])-273.15,c(-1:-3262),col=QG*900)    
# }
# legend("topright",legend=Qlist,col=Qlist*900,lty=1,ncol=2)
# 
# plot(EDCage,type="l",ylim=c(0e5,11e5))
# for(QG in Qlist){
#   age_load<-as.matrix(read.table(sprintf("export/%s/ageRel_m_%.2f_QG_%2.1f_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
#   lines(3248-age_load[,1],age_load[,2],col=QG*1000,type="l")
# }
# legend("topleft",legend=Qlist,col=Qlist*1000,lty=1,ncol=2)
# 
# plot(EDCage,type="l",xlim=c(2300,3300),ylim=c(2e5,11e5))
# for(QG in Qlist){
#   age_load<-as.matrix(read.table(sprintf("export/%s/ageRel_m_%.2f_QG_%2.1f_CN.dat",folder,m,QG*1000),header=FALSE, sep="\t"))
#   lines(3248-age_load[,1],age_load[,2],col=QG*1000,type="l")
# }
# legend("topleft",legend=Qlist,col=Qlist*1000,lty=1,ncol=2)
# lines(EDCage,type="l",lwd=2)
# 
# #V Fischer 0.4, 55-57 K Paterson
# 
# #V Fischer 0.4-0.45, 53-55 K Sc WRONG E
# 
# #V Parenin 3-4, 57-59 K Paterson 
# 
# #V Parenin 2.5-4, 55-57 K shwander WRONG E
# 
# #Without SE:
# 
# #V Parenin 3-3.5, 61-63 K Paterson To rerun
# #Add m=2.5 and m=2
# 
# #V Fischer 0.4-0.5, 58-60 K 
# 
# 
# temp_load<-as.matrix(read.table("export/temp_profile_m_0.20_QG_50_CN.dat",header=FALSE, sep="\t"))
# age_load<-as.matrix(read.table("export/ageRel_m_0.20_QG_50_CN.dat",header=FALSE, sep="\t"))
# 
# 
# 
# temp_load<-as.matrix(read.table("export/temp_profile_m_0.20_QG_50_CN.dat",header=FALSE, sep="\t"))
# temp_load2<-as.matrix(read.table("export/temp_profile_m_0.50_QG_45_CN.dat",header=FALSE, sep="\t"))
# age_load<-as.matrix(read.table("export/ageRel_m_0.20_QG_50_CN.dat",header=FALSE, sep="\t"))
# age_load2<-as.matrix(read.table("export/ageRel_m_0.20_QG_45_CN.dat",header=FALSE, sep="\t"))
# 
# par(mfrow=c(1,1))
# 
# plot(EDCf,type="l")
# lines(rev(temp_load[1:3262])-273.15,c(-1:-3262),col="red")
# lines(rev(temp_load2[1:3262])-273.15,c(-1:-3262),col="blue")
# 
# plot(EDCage,type="l")
# lines(3248-age_load[,1],age_load[,2],col="red",type="b")
# lines(3248-age_load2[,1],age_load2[,2],col="blue",type="b")
# 
# 
# 
# 
# 
# 
# 
# 
# plot(summary1[,1],(summary1[,4]),col=colo[100*summary1[,2]])
# plot(summary1[,2],(summary1[,4]),col=colo[10*summary1[,4]])
# 
# summary1[,c(2,3)]
# 
# ,col=colo[(summary1[,3]-min(summary1[,3]))/(max(summary1[,3])-min(summary1[,3]))*10])
# 
# (summary1[,3]-min(summary1[,3]))/(max(summary1[,3])-min(summary1[,3]))*10
# 
# colo[1.4]
# 
# summary1
# 
# 
# mean(abs(EDCf[,1]-tempCN[c(8:3262),1]))
# 
# 
# 
# tail(tempCN)
# 
# (EDCf[])
# 
# 
# 
# 

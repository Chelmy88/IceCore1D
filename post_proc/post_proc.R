setwd("~/Documents/KUP/IceCore1D/post_proc")
ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}

### tcor2 : 2,2.5,3.5 XX tcor=-1, -1.5 , pcor -5,-10,-15, QG 51-->59, mw 0.45,0.5,0.55
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
running<-function(table){
  len=length(table)
  mel=table
  for (i in 50:(len-50)){
    table[i]=sum(mel[(i-49):(i+49)])/99
  }

  table[1:49]=NA
  table[(len-49):len]=NA
  return (table)
}



read_init<-function(output)
{
  fileName <- paste0("../output/",output,"/init.txt")
  conn <- file(fileName,open="r")
  lines <-readLines(conn)
  init=list()
  for(l in lines)
  {
    l=trim(l)
    if(startsWith(l,"#") || l=="") next
    l=trim(strsplit(gsub("#.*","",l),":")[[1]])
    l=strsplit(l," ")
    init[[l[[1]]]]=l[[2]]
  }
  close(conn)
  return(init)
}

init
read_data<-function(output,init)
{
  data=list()
  data$temperature=data.frame()
  data$age=data.frame()
  for(name in c("QG","MW","PCOR","PCOR","TCOR2","TCOR","DELTAH","FLAT","LEN")){
    data[[name]]=list()
  }
  columns=1
  for (QG in as.numeric(init$QG)){
    for (MW in as.numeric(init$MW)){
      for (PCOR in as.numeric(init$PCOR)){
        for (TCOR in as.numeric(init$TCOR)){
          for (TCOR2 in as.numeric(init$TCOR2)){
            for (DELTAH in as.numeric(init$DELTAH)){
              for (FLAT in as.numeric(init$FLAT)){
                for (LEN in as.numeric(init$LEN)){
                  dirname<-print(sprintf("m_%.3f_Q_%.2f_Pcor_%.0f_Tcor_%.1f_Tcor2_%.1f_dH_%.0f_len_%.0f_flat_%.0f_EDC_Rho_Snow_%s_Thermal_Ice_%s_Thermal_Firn_%s_Heat_Capacity_%s_Rho_Firn_%s_Internal_Energy_%s_Scheme_%s",
                                      MW,QG*1000,PCOR,TCOR,TCOR2,DELTAH,LEN,FLAT,init$RHO_SNOW,init$THERMAL_ICE,init$THERMAL_FIRN,
                                      init$HEAT_CAPACITY,init$RHO_FIRN,init$INTERNAL_ENERGY,init$SCHEME))
                  temperature=as.data.frame(read.table(paste0("../output/",output,"/",dirname,"/temp_profile.dat"),header=FALSE,sep="\t"))
                  temperature[temperature==0]=NA
                  age=as.data.frame(read.table(paste0("../output/",output,"/",dirname,"/age_profile.dat"),header=FALSE,sep="\t"))
                  if(ncol(data$temperature)==0){
                    data$temperature=cbind(c(1:length(temperature[,1])-1),temperature)
                    colnames(data$temperature)=c("depth",dirname)
                    data$age=age
                    colnames(data$temperature)=c("depth",dirname)
                  } else{
                    cnames=c(colnames(data$temperature),dirname)
                    data$temperature=cbind(data$temperature,temperature)
                    colnames(data$temperature)=cnames
                    data$age=cbind(data$age,age[,2])
                    colnames(data$age)=cnames
                  }
                  #Update current column count
                  columns=columns+1
                  #Keep track of column number corresponding to unique values of input parameters
                  if(is.null(data[["QG"]][[sprintf("%.3f",QG)]])){data[["QG"]][[sprintf("%.3f",QG)]]=columns
                  }else{data[["QG"]][[sprintf("%.3f",QG)]]=c(data[["QG"]][[sprintf("%.3f",QG)]],columns)}
                  for(name in c("MW","PCOR","TCOR","TCOR2","DELTAH","FLAT","LEN")){
                    if(is.null(data[[name]][[as.character(eval(parse(text = name)))]])){data[[name]][[as.character(eval(parse(text = name)))]]=columns
                    }else{data[[name]][[as.character(eval(parse(text = name)))]]=c(data[[name]][[as.character(eval(parse(text = name)))]],columns)}
                  }

                }
              }
            }
          }
        }
      }
    }
  }

  return(data)
}




plot_scatter <- function(data,EDC_temp,age_borehole) {
  par(mfrow=c(3,2))
  for(var in c("QG","MW","PCOR","TCOR","TCOR2")){
    col=c()
    mae=c()
    mre=c()
    rmse=c()
    cl=1
    for(ind in sort(names(data[[var]])))
    {
      values=data[[var]][[ind]]
      col=c(col,rep(cl,length(values)))
      for(val in values){
        mae=c(mae,mean(abs(data$temperature[c(1:length(EDC_temp)),val]-EDC_temp)))
        mre=c(mre,mean((abs(data$age[,val]-age_borehole_interpol$y)/
                          age_borehole_interpol$y)[1:200],na.rm=T))
        rmse=c(rmse,sqrt(mean(((data$age[,val]-age_borehole_interpol$y)[1:200])^2,
                              na.rm=T)))
      }
      cl=cl+1
    }
    plot(mae,mre,col=col,main=paste(var,"MRE"))
    legend("top",legend=sort(names(data[[var]])),col=c(1:length(names(data[[var]]))),
           lty=1,bty="n",ncol=length(names(data[[var]])))
  }

}

plot_diff <- function(data,EDC_temp,age_borehole_interpol) {
  par(mfrow=c(5,2))
  for(var in c("QG","MW","PCOR","TCOR","TCOR2")){
    cl=1
    plot(1,xlim=c(-2,2),ylim=c(0,3200),type="n",main=paste(var,"Temperature"))
    for(ind in names(data[[var]]))
    {
      values=data[[var]][[ind]]
      for(val in values){
        lines(data$temperature[c(1:length(EDC_temp)),val]-EDC_temp,1:length(EDC_temp),
              col=cl)
      }
      cl=cl+1
    }
    legend("top",legend=names(data[[var]]),col=c(1:length(names(data[[var]]))),
           lty=1,bty="n",ncol=length(names(data[[var]])))

    cl=1
    plot(1,xlim=c(-0.2,0.2),ylim=c(0,3200),type="n",main=paste(var,"Age"))
    for(ind in names(data[[var]]))
    {
      values=data[[var]][[ind]]
      for(val in values){
        lines((data$age[,val]-age_borehole_interpol$y)/age_borehole_interpol$y,
              age_borehole_interpol$x,col=cl)
      }
      cl=cl+1
    }
    legend("top",legend=names(data[[var]]),col=c(1:length(names(data[[var]]))),
           lty=1,bty="n",ncol=length(names(data[[var]])))
  }
}

EDC_temp<-as.matrix(read.table("../time_series/EDC_temp_forC.dat",header=FALSE,sep="\t"))
age_borehole<-as.matrix(read.table("../time_series/EDC_AICC2012_official.csv",header=TRUE,sep="\t"))
age_borehole<-age_borehole[,c(1,2)]
age_borehole[,1]<-3275-age_borehole[,1]
age_borehole[,2]<-age_borehole[,2]+50

names(data[["MW"]])

output="output_test_CP"
init=read_init(output)
data=read_data(output,init)
age_borehole_interpol=approx(age_borehole[,1],age_borehole[,2],data$age$depth)

pdf(paste0(output,".pdf"), width=8, height=12)
plot_scatter(data,EDC_temp,age_borehole)
plot_diff(data,EDC_temp,age_borehole_interpol)
dev.off()

plot(EDC_temp-data$temperature[1:3276,1])

head(data$temperature[,1])
head(EDC_temp)

colnames(data$temperature)[2]-

plot(EDC_temp,1:length(EDC_temp),type="l",xlim=c(210,230),ylim=c(2800,3300))
lines(data$temperature[,60],data$temperature$depth,col=2)

plot(age_borehole,type="l",xlim=c(3000,3200),ylim=c(0,1e4))
lines(data$age$depth,data$age[,2])

plot(temperature,9.828 * exp(-0.0057 * temperature),type="l")
lines(temperature, 2.22 * (1 - 0.0067 * (temperature - 273.15)))

#
# data$age[,val]-age_borehole_interpol$y
#
#
# data$LEN
#
#
#
# var
#
# mean(abs(data$temperature[c(1:length(EDC_temp)),3]-EDC_temp))
#
# mean(abs(data$age[,3]-age_borehole_interpol$y)/age_borehole_interpol$y,na.rm=T)
#
# sqrt(mean(((data$age[,3]-age_borehole_interpol$y))^2,na.rm=T))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# EDC_temp<-as.matrix(read.table("../time_series/EDC_temp_forC.dat",header=FALSE,sep="\t"))
# age_borehole<-as.matrix(read.table("../time_series/EDC_AICC2012_official.csv",header=TRUE,sep="\t"))
# age_borehole<-age_borehole[,c(1,2)]
# age_borehole[,1]<-3275-age_borehole[,1]
# age_borehole[,2]<-age_borehole[,2]+50
# #age_boreholeC<-as.matrix(read.table("../time_series/EDC_age_forC.dat",header=FALSE,sep="\t"))
#
# surface_temp<-as.matrix(read.table("../time_series/LR04-EDC_temp_1Myr.dat",header=FALSE,sep="\t"))
#
# surface_temp_ma=ma(surface_temp)
# surface_temp_ma[which(is.na(surface_temp_ma))]=surface_temp[which(is.na(surface_temp_ma))]
# surface_temp_ma[10000]=(surface_temp_ma[9999]+surface_temp_ma[10001])/2
# surface_temp_ma_cor=surface_temp_ma
# surface_temp_ma_cor[9880:9998]=surface_temp_ma_cor[9880:9998]-0.3
# plot(surface_temp[9800:10001]-273.15,type="l")
# lines(surface_temp_ma[9800:10001]-273.15,type="l",col=2)
# lines(surface_temp_ma_cor[9800:10001]-273.15,type="l",col=3)
# write.table(surface_temp_ma,"../time_series/LR04-EDC_temp_1Myr_ma.dat",sep="\t",row.names = F,
#             col.names = F)
# write.table(surface_temp_ma_cor,"../time_series/LR04-EDC_temp_1Myr_ma_cor.dat",sep="\t",row.names = F,
#             col.names = F)
#
#
#
#
#
# expand_var <- function(raw_data,var,var_value,conc_res,level)
# {
#   variables=c("QG","MW","PCOR","TCOR2","TCOR","DELTAH","FLAT","LEN")
#   var_level=which(variables==var)
#   var_level
#   if(level==var_level)
#   {
#     raw_data=raw_data[[var_value]]
#   }
#   if(is.null(names(raw_data)))
#   {
#
#     conc_res=cbind(conc_res,raw_data)
#   }
#   else {
#     level=level+1
#     for(n in names(raw_data))
#     {
#       conc_res=expand_var(raw_data[[n]],var,var_value,conc_res,level)
#     }
#   }
#   return(conc_res)
# }
#
#
# summarize_by_var<-function(raw_data,var,range,init)
# {
#   data=list()
#   for (val in init[[var]])
#   {
#     data[[val]]=matrix(nrow=3400, ncol=1)
#     data[[val]]=expand_var(raw_data,var,val,data[[val]],1)[range,]
#     data[[val]]=sweep(data[[val]],1,as.vector(EDC_temp),"-")
#     data[[val]]<-data[[val]][,-1]
#   }
#   return(data)
# }
#
# summarize_by_var_raw<-function(raw_data,var,range,init)
# {
#   data=list()
#   for (val in init[[var]])
#   {
#     data[[val]]=matrix(nrow=3400, ncol=1)
#     data[[val]]=expand_var(raw_data,var,val,data[[val]],1)[range,]
#     data[[val]]<-data[[val]][,-1]
#   }
#   return(data)
# }
# var="QG"
# val=init[[var]]
# plot(data$age$`0.054`$`0.5`$`5`$`2`$`0`$`100`$`500`$`5000`)
#
# range=c(1:3276)
# data$age
#
# summarize_by_var_age<-function(raw_data,var,range,init)
# {
#   data=list()
#   for (val in init[[var]])
#   {
#     data[[val]]=matrix(nrow=680, ncol=1)
#     data[[val]]=expand_var(raw_data,var,val,data[[val]],1)[range,]
#     data[[val]]<-data[[val]][,-1]
#   }
#   return(data)
# }
#
#
#
#
#
#
# data$age
#
# data_var=summarize_by_var_raw(data$data,"QG",1:3276,init)
# data_age=summarize_by_var_age(data$age, "QG",1:3276,init)
#
# plot(data_var[[1]][,1],c(1:length(data_var[[1]][,1])),type="l")
# lines(data_var[[1]][,2],c(1:length(data_var[[1]][,1])),type="l",col=2)
#
#
#
#
#
# output="cor_SC_WE_ADD_350"
#
# init=read_init(output)
#
# init
#
# data=read_data(output,init)
#
# variables=c("QG","MW","PCOR","TCOR2","TCOR","DELTAH","FLAT","LEN")
# colors=c("#0033FF55","#00FF0055","#FF000055","#FF222255","#88888855","#FFFF0055")
#
# par(mar=c(3,3.5,2,0.5),mfrow=c(2,2))
# for (v in variables){
#   data_var=summarize_by_var(data$data,v,1:3276,init)
#
#   plot(EDC_temp-EDC_temp,1:3276,type="l",col="white",xlab="",ylab="",ylim=c(0,3200),xlim=c(-1,1),main="Temperature difference")
#   col=1
#   for (n in names(data_var)){
#     d=data_var[[n]]
#     polygon(c(d[,which(colSums(d,na.rm=T)==min(colSums(d,na.rm=T)))],
#               rev(d[,which(colSums(d,na.rm=T)==max(colSums(d,na.rm=T)))])),
#             c(1:3276,rev(1:3276)),col=colors[col],border=NA)
#     col=col+1
#   }
#   legend("topright",legend = names(data_var),title=v,fill=colors[c(1:col)])
# }
#
#
#
# par(mar=c(3,3.5,2,0.5),mfrow=c(2,2))
# for (v in variables){
#   data_var=summarize_by_var(data,v,1:3276,init)
#   for (n in names(data_var)){
#     d=data_var[[n]]
#     head(d)
#     lines(data$`0.050`$`0.4`$`0`$`0`$`1`$`100`$`5000`$`5000`,1:3400)
#     col=col+1
#   }
#   legend("topright",legend = names(data_var),title=v,fill=colors[c(1:col)])
# }
#
# length(data$`0.050`$`0.4`$`0`$`0`$`1`$`100`$`5000`$`5000`)
#
# for (v in variables){
#   v="QG"
#   data_var=summarize_by_var(data,v,1:3276,init)
#   col=1
#   plot(x=c(0,0.15),y=c(0,30),type="n")
#   df <- data.frame(dummy=c(1:512))
#   for (n in names(data_var)){
#     d=data_var[[n]]
#     den=density(sqrt(colSums(d*d)/length(d)))
#     df_tmp=data.frame(den$x,den$y)
#     colnames(df_tmp)=c(paste0("c_",n,"_x"),paste0("c_",n,"_y"))
#     df=cbind(df,df_tmp)
#     #lines(density(sqrt(colSums(d*d)/length(d))),col=colors[col])
#     #fig <- fig %>% add_trace(density(sqrt(colSums(d*d)/length(d))),y = ~y,x=~x, name = 'trace 0',mode = 'lines')
#     col=col+1
#   }
#   legend("topright",legend = names(data_var),title=v,lty=1,col=colors[c(1:col)])
# }
#
# fig=plot_ly(df,type="scatter")
# for (n in names(data_var)){
#   fig <- fig %>% add_trace(x = df[[paste0("c_",n,"_x")]], y= df[[paste0("c_",n,"_y")]], name = n, mode = 'lines')
# }
#
# fig2=plot_ly(df,type="scatter")
# for (n in names(data_var)){
#   fig2 <- fig2 %>% add_trace(x = df[[paste0("c_",n,"_x")]], y= df[[paste0("c_",n,"_y")]], name = n, mode = 'lines')
# }
# fig <- subplot(fig, fig2)
# fig
#
# names(data_var)
#
#
#
# lines(density(sqrt(colSums(d*d)/length(d))))
# lines(den$x,den$y)
# fig <- plot_ly(x=den$x,y=den$y, name = 'trace 0',mode = 'lines')
# fig <- fig %>% add_trace(x=den$x,y=den$y, name = 'trace 1',mode = 'scatter',col="red")
# fig
#
#
# par(mar=c(3,3.5,2,0.5),mfrow=c(2,2))
# for (v in variables){
#   data_var=summarize_by_var_age(age,v,1:650,init)
#   plot(age_borehole[,c(2,1)],type="l",col=1,xlab="",ylab="",main="Age")
#   col=1
#   for (n in names(data_var)){
#     d=data_var[[n]]
#     polygon(c(d[,which(colSums(d)==min(colSums(d)))],
#               rev(d[,which(colSums(d)==max(colSums(d)))])),
#             c(5*(1:650),5*rev(1:650)),col=colors[col],border=NA)
#     col=col+1
#   }
#   legend("topright",legend = names(data_var),title=v,fill=colors[c(1:col)])
# }
# library(plotly)
# plot_ly(x=c(1:650),y=data_var[[n]][,1],type="scatter")
#
# plot(c(1:650),data_var[[n]][,1])
#
# data_var[[n]][,1]
#
#
# xlength(data_var[[n]][,1])


#
# setThermalIce_CP<-function(temperature) {
#   K = 9.828 * exp(-0.0057 * temperature);
#   return(K)
# }
#
# setThermalIce_GO<-function(temperature) {
#   K = 2.22 * (1 - 0.0067 * (temperature - 273.15));
#   return(K)
# }
#
# setThermalFirn_CP<-function(K,rho,cp) {
#   K = 2. * K * rho / (3 * 923 - rho)
#   return(K)
# }
#
# setThermalFirn_SC<-function(K, rho,cp) {
#   K = K * (rho / 923)**(2 - 0.5 * rho / 923)
#   return(K)
# }
#
# setThermalFirn_CP_LIN<-function(K,rho,cp) {
#   K=setThermalFirn_CP(K,rho,cp)
#   k0 = 25 * density[3275] * cp[3275] / (365.25 * 24 * 3600);
#   k100 = K[3275 - 100];
#   for (li in c(3275:(3275 - 100))) {
#     K[li] = k0 + (3275 - li) / 100. * (k100 - k0);
#   }
#   return(K)
# }
#
# setThermalFirn_SC_LIN<-function(K,rho,cp) {
#   K=setThermalFirn_SC(K,rho,cp)
#   k0 = 25 * density[3275] * cp[3275] / (365.25 * 24 * 3600);
#   k100 = K[3275 - 100];
#   for (li in c(3275:(3275 - 100))) {
#     K[li] = k0 + (3275 - li) / 100. * (k100 - k0);
#   }
#   return(K)
# }
#
# setThermalFirn_CP_AL<-function(K,rho,rho_ice,temp,cp) {
#   K1=setThermalFirn_CP(K,rho,cp)
#   K2=rho*cp*(1 - 0.00882*(temp + 30-273.15))*(-1.229e-14 * rho * rho * rho +
#                                                 2.1312*1e-11*rho*rho - 9.4e-9*rho +1.779e-6);
#   K=(K2-K1)* ((rho-rho_ice)/(rho[3275]-rho_ice))^1 + K1
#   return(K)
# }
#
# setThermalFirn_SC_AL<-function(K,rho,rho_ice,temp,cp) {
#   K1=setThermalFirn_SC(K,rho,cp)
#   K2=rho*cp*(1 - 0.00882*(temp + 30-273.15))*(-1.229e-14 * rho * rho * rho +
#                                                 2.1312*1e-11*rho*rho - 9.4e-9*rho +1.779e-6);
#   K=(K2-K1)* ((rho-rho_ice)/(rho[3275]-rho_ice))^1 + K1
#   return(K)
# }
#
#
# setThermalFirn_SC_WE_LIN<-function(K,rho,rho_ice,cp) {
#   K1=setThermalFirn_SC(K,rho,cp)
#   K_eff_WE=c(0.32,0.44,0.51,0.64,0.76,0.78,0.77,0.71,0.65,0.68,0.71)
#   K=(K_eff_WE[11]-K1)* ((rho-rho_ice)/(rho[3275-10]-rho_ice))^4 + K1
#   K[3275:3265]=K_eff_WE
#   return(K)
# }
#
# setThermalFirn_CP_WE_LIN<-function(K,rho,rho_ice,cp) {
#   K1=setThermalFirn_CP(K,rho,cp)
#   K_eff_WE=c(0.32,0.44,0.51,0.64,0.76,0.78,0.77,0.71,0.65,0.68,0.71)
#   K=(K_eff_WE[11]-K1)* ((rho-rho_ice)/(rho[3275-10]-rho_ice))^4 + K1
#   K[3275:3265]=K_eff_WE
#   return(K)
# }
#
# setThermalFirn_SC_WE_ADD<-function(K,rho,rho_ice,cp) {
#   K=setThermalFirn_SC(K,rho,cp)
#   K_eff_WE_delta=c(0.056, 0.023, 0.071, 0.154, 0.238, 0.250, 0.196, 0.112, 0.027, 0.014, 0.000)
#   K[3275:3265]=K[3275:3265]+K_eff_WE_delta
#   return(K)
# }
#
# setThermalFirn_CP_WE_ADD<-function(K,rho,rho_ice,cp) {
#   K=setThermalFirn_CP(K,rho,cp)
#   K_eff_WE_delta=c(0.056, 0.023, 0.071, 0.154, 0.238, 0.250, 0.196, 0.112, 0.027, 0.014, 0.000)
#   K[3275:3265]=K[3275:3265]+K_eff_WE_delta
#   return(K)
# }
#
#
#
# setThermalFirn_ST<-function(K, rho,cp) {
#   K = 0.138-1.01*rho/1000+3.233*rho^2/1000/1000
#   return(K)
# }
#
# rho=density
# rho_ice=ice_density
#
# setThermalFirn_CP_ST<-function(K, rho,rho_ice,cp) {
#   K1 = setThermalFirn_CP(K_ice,density,cp)
#   K2 = 0.138-1.01*rho/1000+3.233*rho^2/1000/1000
#   k = rep(0,3275)
#   li=3275
#   while(rho[li]<600)
#   {
#     K[li]=K2[li]
#     li=li-1
#   }
#   for(l in c(li:1))
#   {
#     K[l]=(K2[l]-K1[l])* ((rho[l]-rho_ice[l])/(rho[li]-rho_ice[l]))^1 + K1[l]
#   }
#   return(K)
# }
# setThermalFirn_SC_ST<-function(K, rho,rho_ice,cp) {
#   K1 = setThermalFirn_SC(K_ice,density,cp)
#   K2 = 0.138-1.01*rho/1000+3.233*rho^2/1000/1000
#   k = rep(0,3275)
#   li=3275
#   while(rho[li]<600)
#   {
#     K[li]=K2[li]
#     li=li-1
#   }
#   for(l in c(li:1))
#   {
#     K[l]=(K2[l]-K1[l])* ((rho[l]-rho_ice[l])/(rho[li]-rho_ice[l]))^1 + K1[l]
#   }
#   return(K)
# }
#
#
# setHeatCapacity_CP<-function(temperature){
#   cp = 152.5 + 7.122 * temperature
#   return(cp)
# }
#
# setHeatCapacity_AL<-function(temperature,density,ice_density){
#   cp_sc=152.5 + 7.122 * temperature
#   cp=-27.796+7.7752*temperature
#   cp=(cp-cp_sc)* ((density-ice_density)/(density[3275]-ice_density)) + cp_sc
#   return(cp)
# }
#
#
#
#
#
#
#
# name
#
#
#
# n <- 6
# palette <- distinctColorPalette(n)
#
# THERMAL_ICE=c("CP","GO")
# THERMAL_FIRN=c("CP","SC","CP_WE_ADD","CP_WE_LIN","SC_WE_ADD", "SC_WE_LIN")
# HEAT_CAPACITY=c("CP","CP_AL")
# RHO_SNOW=c("300","350","400")
#
# RHO_SNOW="400"
# data=list()
#
#
#
# name_summary=paste0("Comparison_RHO_SNOW_",RHO_SNOW)
#
# for(TF in THERMAL_FIRN){
#   name=paste0(TF,"_",RHO_SNOW)
#   output=paste0("output/",name,"/")
#   init=read_init(output)
#
#   temp=as.data.frame(read.table(paste0("../",output,"m_",format(as.numeric(init$MW),nsmall=3),"_Q_",format(as.numeric(init$QG)*1000,nsmall=2),"_Pcor_",init$PCOR,"_Tcor_",format(as.numeric(init$TCOR),nsmall=1),"_Tcor2_",format(as.numeric(init$TCOR2),nsmall=1),"_dH_",init$DELTAH,"_len_",init$LEN,"_flat_",init$FLAT,"_EDC_Rho_Snow_",init$RHO_SNOW,"_Thermal_Ice_",init$THERMAL_ICE,"_Thermal_Firn_",init$THERMAL_FIRN,"_Heat_Capacity_",init$HEAT_CAPACITY,"_Rho_Firn_",init$RHO_FIRN,"_Internal_Energy_OFF_Scheme_CN/temp_profile.dat")))
#   temp=temp$V1
#
#   density=as.data.frame(read.table(paste0("../",output,"m_",format(as.numeric(init$MW),nsmall=3),"_Q_",format(as.numeric(init$QG)*1000,nsmall=2),"_Pcor_",init$PCOR,"_Tcor_",format(as.numeric(init$TCOR),nsmall=1),"_Tcor2_",format(as.numeric(init$TCOR2),nsmall=1),"_dH_",init$DELTAH,"_len_",init$LEN,"_flat_",init$FLAT,"_EDC_Rho_Snow_",init$RHO_SNOW,"_Thermal_Ice_",init$THERMAL_ICE,"_Thermal_Firn_",init$THERMAL_FIRN,"_Heat_Capacity_",init$HEAT_CAPACITY,"_Rho_Firn_",init$RHO_FIRN,"_Internal_Energy_OFF_Scheme_CN/density_profile.dat")))
#   density=density$V1
#
#   ice_density=as.data.frame(read.table(paste0("../",output,"m_",format(as.numeric(init$MW),nsmall=3),"_Q_",format(as.numeric(init$QG)*1000,nsmall=2),"_Pcor_",init$PCOR,"_Tcor_",format(as.numeric(init$TCOR),nsmall=1),"_Tcor2_",format(as.numeric(init$TCOR2),nsmall=1),"_dH_",init$DELTAH,"_len_",init$LEN,"_flat_",init$FLAT,"_EDC_Rho_Snow_",init$RHO_SNOW,"_Thermal_Ice_",init$THERMAL_ICE,"_Thermal_Firn_",init$THERMAL_FIRN,"_Heat_Capacity_",init$HEAT_CAPACITY,"_Rho_Firn_",init$RHO_FIRN,"_Internal_Energy_OFF_Scheme_CN/pure_ice_density_profile.dat")))
#   ice_density=ice_density$V1
#
#   age=as.data.frame(read.table(paste0("../",output,"m_",format(as.numeric(init$MW),nsmall=3),"_Q_",format(as.numeric(init$QG)*1000,nsmall=2),"_Pcor_",init$PCOR,"_Tcor_",format(as.numeric(init$TCOR),nsmall=1),"_Tcor2_",format(as.numeric(init$TCOR2),nsmall=1),"_dH_",init$DELTAH,"_len_",init$LEN,"_flat_",init$FLAT,"_EDC_Rho_Snow_",init$RHO_SNOW,"_Thermal_Ice_",init$THERMAL_ICE,"_Thermal_Firn_",init$THERMAL_FIRN,"_Heat_Capacity_",init$HEAT_CAPACITY,"_Rho_Firn_",init$RHO_FIRN,"_Internal_Energy_OFF_Scheme_CN/age_profile.dat")))
#
#   if(init$HEAT_CAPACITY=="CP")
#   {
#     cp=setHeatCapacity_CP(temp)
#   } else if(init$HEAT_CAPACITY=="AL")
#   {
#     cp=setHeatCapacity_AL(temp,density,ice_density)
#   }
#   print(init$THERMAL_ICE)
#
#   if(init$THERMAL_ICE=="CP")
#   {
#     print("Using themal ice CP")
#     K_ice=setThermalIce_CP(temp)
#   } else if(init$THERMAL_ICE=="GO")
#   {
#     K_ice=setThermalIce_GO(temp)
#   }
#
#   if(init$THERMAL_FIRN=="CP")
#   {
#     K=setThermalFirn_CP(K_ice,density,cp)
#   } else if(init$THERMAL_FIRN=="SC")
#   {
#     K=setThermalFirn_SC(K_ice,density,cp)
#   } else if(init$THERMAL_FIRN=="CP_LIN")
#   {
#     K=setThermalFirn_CP_LIN(K_ice,density,cp)
#   } else if(init$THERMAL_FIRN=="SC_LIN")
#   {
#     K=setThermalFirn_SC_LIN(K_ice,density,cp)
#   } else if(init$THERMAL_FIRN=="CP_AL")
#   {
#     K=setThermalFirn_CP_AL(K_ice,density,ice_density,temp,cp)
#   } else if(init$THERMAL_FIRN=="SC_AL")
#   {
#     K=setThermalFirn_SC_AL(K_ice,density,ice_density,temp,cp)
#   } else if(init$THERMAL_FIRN=="CP_ST")
#   {
#     K=setThermalFirn_CP_ST(K_ice,density,ice_density,cp)
#   } else if(init$THERMAL_FIRN=="SC_ST")
#   {
#     K=setThermalFirn_SC_ST(K_ice,density,ice_density,cp)
#   }
#   else if(init$THERMAL_FIRN=="CP_WE_ADD")
#   {
#     K=setThermalFirn_CP_WE_ADD(K_ice,density,ice_density,cp)
#   }
#   else if(init$THERMAL_FIRN=="CP_WE_LIN")
#   {
#     K=setThermalFirn_CP_WE_LIN(K_ice,density,ice_density,cp)
#   }
#   else if(init$THERMAL_FIRN=="SC_WE_ADD")
#   {
#     K=setThermalFirn_SC_WE_ADD(K_ice,density,ice_density,cp)
#   }
#   else if(init$THERMAL_FIRN=="SC_WE_LIN")
#   {
#     K=setThermalFirn_SC_WE_LIN(K_ice,density,ice_density,cp)
#   }
#
#   heat_diff=(K/density/cp)*365.25*24*3600
#
#   data[[init$RHO_SNOW]][[init$THERMAL_ICE]][[init$THERMAL_FIRN]][[init$HEAT_CAPACITY]]=list("temperature"=temp,"rho"=density,"K_ice"=K_ice,"K"=K,"CP"=cp,"heat_diff"=heat_diff,"age"=age)
#
#   pdf(paste0("../output/",name,".pdf"),width=12,height=8)
#   par(mfrow=c(3,4),mar=c(3,3,2,1),mgp=c(1.7,0.5,0),oma=c(0,0,2,0))
#
#   plot(EDC_temp,1:3276,type="l",xlab="T (K)",ylab="Height (m)",ylim=c(2800,3300),xlim=c(217,224),main="Temperature Top",col=2)
#   lines(243.125+c(3150:3300)* -0.007575,c(3150:3300),lty=2)
#   lines(temp[c(2800:3275)],c(2800:3275))
#
#   plot(EDC_temp,1:3276,type="l",ylim=c(0,3300),main="Temperature",col=2,ylab="Height (m)",xlab="T (K)",xlim=c(215,275))
#   lines(temp[c(1:3275)],c(1:3275),ylab="Height (m)")
#
#   plot(temp[c(1:3276)]-EDC_temp,1:3276,type="l",ylim=c(0,3300),main="Temperature difference",col=2,ylab="Height (m)",xlab="Delta T (K)",xlim=c(-1,1))
#   lines(temp[c(1:3276)]-c(EDC_temp[1:3149],243.125+c(3150:3276)* -0.007575),c(1:3276),ylab="Height (m)")
#
#   plot(temp[c(1:3276)]-EDC_temp,1:3276,type="l",ylim=c(2800,3300),main="Temperature difference top",col=2,ylab="Height (m)",xlab="Delta T (K)",xlim=c(-1,1))
#   lines(temp[c(1:3276)]-c(EDC_temp[1:3149],243.125+c(3150:3276)* -0.007575),c(1:3276),ylab="Height (m)")
#
#
#   plot(cp[c(2800:3275)],c(2800:3275),type="l",main="Heat capacity",ylab="Height (m)",xlab="CP (J/kg/K)",xlim=c(1700,1750))
#   plot(cp[c(1:3275)],c(1:3275),type="l",main="Heat capavity top",ylab="Height (m)",xlab="CP (J/kg/K)",xlim=c(1700,2100))
#
#
#   plot(K_ice[c(2800:3275)],c(2800:3275),type="l",main="Pure ice heat conductivity top",ylab="Height (m)",xlab="K (W/m/K)",xlim=c(2.75,2.83))
#   plot(K_ice[c(1:3275)],c(1:3275),type="l",main="Pure ice Heat conductivity",ylab="Height (m)",xlab="K (W/m/K)",xlim=c(2,2.9))
#
#
#   plot(K[c(2800:3275)],c(2800:3275),type="l",main="Heat conductivity top",ylab="Height (m)",xlab="K (W/m/K)",xlim=c(0,3.5))
#   plot(K[c(1:3275)],c(1:3275),type="l",main="Heat conductivity top",ylab="Height (m)",xlab="K (W/m/K)",xlim=c(0,3.5))
#
#   plot(heat_diff[2800:3275],c(2800:3275),type="l",main="Heat diffusivity top",ylab="Height (m)",xlab="a (m2/yr)",xlim=c(5,60))
#   plot(heat_diff[1:3275],c(1:3275),type="l",main="Heat diffusivity",ylab="Height (m)",xlab="a (m2/yr)",xlim=c(5,60))
#
#   title(paste("Rho snow:",init$RHO_SNOW,"Thermal ice:",init$THERMAL_ICE,"Thermal firn:",init$THERMAL_FIRN,"Heat capacity:",init$HEAT_CAPACITY),outer=T)
#   dev.off()
# }
#
#
# pdf(paste0("../output/",name_summary,".pdf"),width=12,height=8)
#
# library(randomcoloR)
#
# params=names(data[[RHO_SNOW]]$CP)
#
# par(mfrow=c(3,4),mar=c(3,3,2,1),mgp=c(1.7,0.5,0),oma=c(0,0,2,0))
#
# plot(EDC_temp,1:3276,type="l",xlab="T (K)",ylab="Height (m)",ylim=c(3175,3300),xlim=c(217,224),main="Temperature Top",col=1)
# lines(243.125+c(3150:3300)* -0.007575,c(3150:3300),lty=2)
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$temperature[c(2800:3275)],c(2800:3275),col=palette[i])
#   i=i+1
# }
# legend("topright",legend=c("measures",params),col=c(1,palette),lty=1,bty="n",cex=0.9)
#
#
# plot(EDC_temp,1:3276,type="l",ylim=c(0,3300),main="Temperature",col=2,ylab="Height (m)",xlab="T (K)")
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$temperature[c(1:3275)],c(1:3275),col=palette[i])
#   i=i+1
# }
# legend("topright",legend=c("measures",params),col=c(1,palette),lty=1,bty="n",cex=0.9)
#
#
# plot(temp[c(1:3276)]-EDC_temp,1:3276,type="n",ylim=c(0,3300),main="Temperature difference",col=2,ylab="Height (m)",xlab="Delta T (K)",xlim=c(-0.5,1))
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$temperature[c(1:3276)]-c(EDC_temp[1:3149],243.125+c(3150:3276)* -0.007575),c(1:3276),ylab="Height (m)",col=palette[i])
#
#   i=i+1
# }
# legend("topleft",legend=c(params),col=palette,lty=1,bty="n",cex=0.9)
#
#
# plot(temp[c(3175:3276)]-EDC_temp[c(3175:3276)],3175:3276,type="n",ylim=c(3175,3300),main="Temperature difference top",col=2,ylab="Height (m)",xlab="Delta T (K)",xlim=c(-0.5,1))
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$temperature[c(1:3276)]-c(EDC_temp[1:3149],243.125+c(3150:3276)* -0.007575),c(1:3276),ylab="Height (m)",col=palette[i])
#
#   i=i+1
# }
# legend("topleft",legend=c(params),col=palette,lty=1,bty="n",cex=0.9)
#
#
# plot(K[c(3175:3275)],c(3175:3275),type="n",main="Heat conductivity top",ylab="Height (m)",xlab="K (W/m/K)",xlim=c(0,3))
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$K[c(1:3275)],c(1:3275),col=palette[i])
#   i=i+1
# }
# legend("bottomleft",legend=c(params),col=palette,lty=1,bty="n",cex=0.9)
#
#
# plot(K[c(1:3275)],c(1:3275),type="l",main="Heat conductivity",ylab="Height (m)",xlab="K (W/m/K)",xlim=c(0,3))
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$K[c(1:3275)],c(1:3275),col=palette[i])
#   i=i+1
# }
# legend("bottomleft",legend=c(params),col=palette,lty=1,bty="n",cex=0.9)
#
#
# plot((heat_diff)[3175:3275],c(3175:3275),type="n",main="Heat diffusivity top",ylab="Height (m)",xlab="a (m2/yr)",xlim=c(5,60))
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$heat_diff[c(1:3275)],c(1:3275),col=palette[i])
#   i=i+1
# }
# legend("bottomleft",legend=c(params),col=palette,lty=1,bty="n",cex=0.9)
#
#
# plot((heat_diff)[1:3275],c(1:3275),type="n",main="Heat diffusivity",ylab="Height (m)",xlab="a (m2/yr)",xlim=c(5,60))
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$heat_diff[c(1:3275)],c(1:3275),col=palette[i])
#   i=i+1
# }
# legend("bottomleft",legend=c(params),col=palette,lty=1,bty="n",cex=0.9)
#
# plot(age_borehole[,2],age_borehole[,1],type="l",ylab="Hieght (m)",xlab="Age (yr)",main="Ice age")
# i=1
# for(p in params){
#   lines(data[[RHO_SNOW]]$CP[[p]]$CP$age$V2,data[[RHO_SNOW]]$CP[[p]]$CP$age$V1,col=palette[i])
#   i=i+1
# }
# legend("bottomleft",legend=c(params),col=palette,lty=1,bty="n",cex=0.9)
# dev.off()
#
#
#
#
#
# ########################################
# pdf(paste0("../output/parameters.pdf"),width=12,height=8)
#
#
# par(mfrow=c(2,2),mar=c(3,3,2,1),mgp=c(1.7,0.5,0),oma=c(0,0,2,0))
#
# rho=density
# cp = 152.5 + 7.122 * temp
#
# # cp2 = 2096 + (7.7752*(temp-273.15))
# # cp3=setHeatCapacity_AL(temp,density,ice_density)
# #
# # plot(cp[3000:3275],c(3000:3275),type="l",xlim=c(1650,1750),ylab="Depth",xlab="CP", main="Heat capacity comparison top")
# # lines(cp2[3000:3275],c(3000:3275),col=2)
# # lines(cp3[3000:3275],c(3000:3275),col=3)
# # legend("topright",legend=c("CP SC","CP AL","CP AL to CP"),col=c(1,2,3),lty=1,bty="n")
# #
# # plot(cp[1:3275],c(1:3275),type="l",xlim=c(1600,2100),ylab="Depth",xlab="CP", main="Heat capacity comparison")
# # lines(cp2[1:3275],c(1:3275),col=2)
# # lines(cp3[1:3275],c(1:3275),col=3)
# # legend("topright",legend=c("CP SC","CP AL","CP AL to CP"),col=c(1,2,3),lty=1,bty="n")
#
#
# # K_ice_GO = setThermalIce_GO(temp)
# # K_ice_CP = setThermalIce_CP(temp)
# #
# # plot(K_ice_CP[3000:3275],c(3000:3275),type="l",xlim=c(2.6,3.3),ylab="Depth",xlab="K", main="K pure ice comparison top")
# # lines(K_ice_GO[1:3275],c(1:3275),col=2)
# # legend("topright",legend=c("K ice CP","K ice GO"),col=c(1,2),lty=1,bty="n")
# #
# # plot(K_ice_CP[1:3275],c(1:3275),type="l",xlim=c(2,3.3),ylab="Depth",xlab="K", main="K pure ice comparison ")
# # lines(K_ice_GO[1:3275],c(1:3275),col=2)
# # legend("topright",legend=c("K ice CP","K ice GO"),col=c(1,2),lty=1,bty="n")
#
# K_SC = setThermalFirn_SC(K_ice,density,cp)
#
# K_CP = setThermalFirn_CP(K_ice,density,cp)
#
# K_SC_LIN = setThermalFirn_SC_LIN(K_ice,density,cp)
#
# K_CP_LIN = setThermalFirn_CP_LIN(K_ice,density,cp)
#
# K_SC_AL =  setThermalFirn_SC_AL(K_ice,density,ice_density,temp,cp)
#
# K_CP_AL =  setThermalFirn_CP_AL(K_ice,density,ice_density,temp,cp)
#
# K_ST =  setThermalFirn_ST(K_ice,density,cp)
#
# K_CP_ST =  setThermalFirn_CP_ST(K_ice,density,ice_density,cp)
#
# K_SC_ST =  setThermalFirn_SC_ST(K_ice,density,ice_density,cp)
#
# K_SC_WE_LIN=setThermalFirn_SC_WE_LIN(K_ice,density,ice_density,cp)
#
# K_CP_WE_LIN=setThermalFirn_CP_WE_LIN(K_ice,density,ice_density,cp)
#
# K_AL=rho*cp*(1 - 0.00882*(temp + 30-273.15))*(-1.229e-14 * density * density * density +
#                                                 2.1312*1e-11*density*density - 9.4e-9*density +1.779e-6)
#
# plot(K_SC[3175:3275],c(3175:3275),type="l",xlim=c(0.3,3.5),ylab="Depth",xlab="K", main="K comparison top")
# lines(K_CP[1:3275],c(1:3275),col=2)
# lines(K_AL[1:3275],c(1:3275),col=3)
# lines(K_ST[1:3275],c(1:3275),col=4,lty=1)
#
# lines(K_SC_LIN[1:3275],c(1:3275),col=1,lty=2)
# lines(K_CP_LIN[1:3275],c(1:3275),col=2,lty=2)
#
# lines(K_SC_AL[1:3275],c(1:3275),col=1,lty=3)
# lines(K_CP_AL[1:3275],c(1:3275),col=2,lty=3)
#
# lines(K_CP_ST[1:3275],c(1:3275),col=1,lty=4)
# lines(K_SC_ST[1:3275],c(1:3275),col=2,lty=4)
#
# lines(K_SC_WE_LIN[1:3275],c(1:3275),col=1,lty=2,lwd=3)
# lines(K_CP_WE_LIN[1:3275],c(1:3275),col=2,lty=2,lwd=3)
#
#
# legend("topright",legend=c("K SC","K CP","K AL","K ST","TO LIN","TO AL","TO ST"),col=c(1,2,3,4,1,1,1),lty=c(1,1,1,1,2,3,4),bty="n")
#
#
#
#
#
#
#
# plot(K_SC[1:3275],c(1:3275),type="l",xlim=c(0.4,3.5),ylab="Depth",xlab="K (W/m/K)", main="K comparison")
# lines(K_CP[1:3275],c(1:3275),col=2)
# lines(K_AL[1:3275],c(1:3275),col=3)
# lines(K_ST[1:3275],c(1:3275),col=4,lty=1)
#
# lines(K_SC_LIN[1:3275],c(1:3275),col=1,lty=2)
# lines(K_CP_LIN[1:3275],c(1:3275),col=2,lty=2)
#
# lines(K_SC_AL[1:3275],c(1:3275),col=1,lty=3)
# lines(K_CP_AL[1:3275],c(1:3275),col=2,lty=3)
#
# lines(K_CP_ST[1:3275],c(1:3275),col=1,lty=4)
# lines(K_SC_ST[1:3275],c(1:3275),col=2,lty=4)
#
# legend("topright",legend=c("K SC","K CP","K AL","K ST","TO LIN","TO AL","TO ST"),col=c(1,2,3,4,1,1,1),lty=c(1,1,1,1,2,3,4),bty="n")
#
#
# heat_diff_SC=(K_SC/density/cp)*365.25*24*3600
# heat_diff_CP=(K_CP/density/cp)*365.25*24*3600
# heat_diff_SC_LIN=(K_SC_LIN/density/cp)*365.25*24*3600
# heat_diff_CP_LIN=(K_CP_LIN/density/cp)*365.25*24*3600
# heat_diff_SC_AL=(K_SC_AL/density/cp)*365.25*24*3600
# heat_diff_CP_AL=(K_CP_AL/density/cp)*365.25*24*3600
# heat_diff_SC_ST=(K_SC_ST/density/cp)*365.25*24*3600
# heat_diff_CP_ST=(K_CP_ST/density/cp)*365.25*24*3600
# heat_diff_AL=(K_AL/density/cp)*365.25*24*3600
# heat_diff_ST=(K_ST/density/cp)*365.25*24*3600
#
#
# plot(heat_diff_SC[3000:3275],c(3000:3275),type="l",ylab="Depth",xlab="Heat diffusivity (m^2/yr)", main="Heat diffusivity top",xlim=c(10,70))
# lines(heat_diff_CP[3000:3275],c(3000:3275),type="l",col=2)
# lines(heat_diff_AL[3000:3275],c(3000:3275),type="l",col=3)
# lines(heat_diff_ST[3000:3275],c(3000:3275),type="l",col=4)
# lines(heat_diff_SC_LIN[3000:3275],c(3000:3275),type="l",col=1,lty=2)
# lines(heat_diff_CP_LIN[3000:3275],c(3000:3275),type="l",col=2,lty=2)
# lines(heat_diff_SC_AL[3000:3275],c(3000:3275),type="l",col=1,lty=3)
# lines(heat_diff_CP_AL[3000:3275],c(3000:3275),type="l",col=2,lty=3)
# lines(heat_diff_SC_ST[3000:3275],c(3000:3275),type="l",col=1,lty=4)
# lines(heat_diff_CP_ST[3000:3275],c(3000:3275),type="l",col=2,lty=4)
# legend("topright",legend=c("K SC","K CP","K AL","K ST","TO LIN","TO AL","TO ST"),col=c(1,2,3,4,1,1,1),lty=c(1,1,1,1,2,3,4),bty="n")
#
# plot(heat_diff_SC[1:3275],c(1:3275),type="l",ylab="Depth",xlab="Heat diffusivity (m^2/yr)", main="Heat diffusivity",xlim=c(10,70))
# lines(heat_diff_CP[1:3275],c(1:3275),type="l",col=2)
# lines(heat_diff_AL[1:3275],c(1:3275),type="l",col=3)
# lines(heat_diff_ST[1:3275],c(1:3275),type="l",col=4)
# lines(heat_diff_SC_LIN[1:3275],c(1:3275),type="l",col=1,lty=2)
# lines(heat_diff_CP_LIN[1:3275],c(1:3275),type="l",col=2,lty=2)
# lines(heat_diff_SC_AL[1:3275],c(1:3275),type="l",col=1,lty=3)
# lines(heat_diff_CP_AL[1:3275],c(1:3275),type="l",col=2,lty=3)
# lines(heat_diff_SC_ST[1:3275],c(1:3275),type="l",col=1,lty=4)
# lines(heat_diff_CP_ST[1:3275],c(1:3275),type="l",col=2,lty=4)
# legend("topright",legend=c("K SC","K CP","K AL","K ST","TO LIN","TO AL","TO ST"),col=c(1,2,3,4,1,1,1),lty=c(1,1,1,1,2,3,4),bty="n")
#
#
#
# dev.off()
#
#
#


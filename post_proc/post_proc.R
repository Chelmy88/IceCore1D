setwd("~/Documents/KUP/IceCore1D/output")

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
EDC_temp<-as.matrix(read.table("../time_series/EDC_temp_forC.dat",header=FALSE,sep="\t"))
age_borehole<-as.matrix(read.table("../time_series/EDC_AICC2012_official.csv",header=TRUE,sep="\t"))
age_borehole<-age_borehole[,c(1,2)]
age_borehole[,1]<-3275-age_borehole[,1]
age_borehole[,2]<-age_borehole[,2]+50
age_boreholeC<-as.matrix(read.table("../time_series/EDC_age_forC.dat",header=FALSE,sep="\t"))

plot(EDC_temp,type="l")

tempeature=list()
age=list()
data=list()

fileName <- "init.txt"
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

for (QG in as.numeric(init$QG)){
  data[[sprintf("%.3f",QG)]]=list()
  for (mw in as.numeric(init$MW)){
    data[[sprintf("%.3f",QG)]][[as.character(mw)]]=list()
    for (pCor in as.numeric(init$PCOR)){
      data[[sprintf("%.3f",QG)]][[as.character(mw)]][[as.character(pCor)]]=list()
      for (tCor2 in as.numeric(init$TCOR2)){
        data[[sprintf("%.3f",QG)]][[as.character(mw)]][[as.character(pCor)]][[as.character(tCor2)]]=list()
        for (tCor in as.numeric(init$TCOR)){
          data[[sprintf("%.3f",QG)]][[as.character(mw)]][[as.character(pCor)]][[as.character(tCor2)]][[as.character(tCor)]]=list()
          for (deltaH in as.numeric(init$DELTAH)){
            data[[sprintf("%.3f",QG)]][[as.character(mw)]][[as.character(pCor)]][[as.character(tCor2)]][[as.character(tCor)]][[as.character(deltaH)]]=list()
            for (flat in as.numeric(init$FLAT)){
              data[[sprintf("%.3f",QG)]][[as.character(mw)]][[as.character(pCor)]][[as.character(tCor2)]][[as.character(tCor)]][[as.character(deltaH)]][[as.character(flat)]]=list()
              for (len in as.numeric(init$LEN)){
                name<-print(sprintf("m_%.3f_Q_%.2f_Pcor_%.0f_Tcor_%.1f_Tcor2_%.1f_dH_%.0f_len_%.0f_flat_%.0f_EDC_Thermal_%s_Firn_%s_Internal_Energy_%s_Scheme_%s",
                                    mw,QG*1000,pCor,tCor,tCor2,deltaH,len,flat,"CP","SC","OFF","CN"))
                data[[sprintf("%.3f",QG)]][[as.character(mw)]][[as.character(pCor)]][[as.character(tCor2)]][[as.character(tCor)]][[as.character(deltaH)]][[as.character(flat)]][[as.character(len)]]<-
                  as.matrix(read.table(paste0(name,"/temperature_today.dat"),header=FALSE,sep="\t"))
                age[[sprintf("%.3f",QG)]][[as.character(mw)]][[as.character(pCor)]][[as.character(tCor2)]][[as.character(tCor)]][[as.character(deltaH)]][[as.character(flat)]][[as.character(len)]]<-
                  as.matrix(read.table(paste0(name,"/age_profile.dat"),header=FALSE,sep="\t"))[,2]
              }
            }
          }
        }
      }
    }
  }
}


expand_var <- function(raw_data,var,var_value,conc_res,level)
{
  variables=c("QG","MW","PCOR","TCOR2","TCOR","DELTAH","FLAT","LEN")
  var_level=which(variables==var)
  if(level==var_level)
  {
    raw_data=raw_data[[var_value]]
  }
  if(is.null(names(raw_data)))
  {

    conc_res=cbind(conc_res,raw_data)
  }
  else {
    level=level+1
    for(n in names(raw_data))
    {
      conc_res=expand_var(raw_data[[n]],var,var_value,conc_res,level)
    }
  }
  return(conc_res)
}

summarize_by_var<-function(raw_data,var,range,init)
{
  data=list()
  for (val in init[[var]])
  {
    data[[val]]=matrix(nrow=3400, ncol=1)
    data[[val]]=expand_var(raw_data,var,val,data[[val]],1)[range,]
    data[[val]]=sweep(data[[val]],1,as.vector(EDC_temp),"-")
    data[[val]]<-data[[val]][,-1]
  }
  return(data)
}

summarize_by_var_age<-function(raw_data,var,range,init)
{
  data=list()
  for (val in init[[var]])
  {
    data[[val]]=matrix(nrow=680, ncol=1)
    data[[val]]=expand_var(raw_data,var,val,data[[val]],1)[range,]
    data[[val]]<-data[[val]][,-1]
  }
  return(data)
}



variables=c("QG","MW","PCOR","TCOR2","TCOR","DELTAH","FLAT","LEN")
colors=c("#0033FF55","#00FF0055","#FF000055","#FF222255","#88888855","#FFFF0055")

par(mar=c(3,3.5,2,0.5),mfrow=c(2,2))
for (v in variables){
  data_var=summarize_by_var(data,v,1:3276,init)
  plot(EDC_temp-EDC_temp,1:3276,type="l",col="white",xlab="",ylab="",ylim=c(0,3200),xlim=c(-1,1),main="Temperature difference")
  col=1
  for (n in names(data_var)){
    d=data_var[[n]]
    polygon(c(d[,which(colSums(d)==min(colSums(d)))],
              rev(d[,which(colSums(d)==max(colSums(d)))])),
            c(1:3276,rev(1:3276)),col=colors[col],border=NA)
    col=col+1
  }
  legend("topright",legend = names(data_var),title=v,fill=colors[c(1:col)])
}

for (v in variables){
  v="QG"
  data_var=summarize_by_var(data,v,1:3276,init)
  col=1
  plot(x=c(0,0.15),y=c(0,30),type="n")
  df <- data.frame(dummy=c(1:512))
  for (n in names(data_var)){
    d=data_var[[n]]
    den=density(sqrt(colSums(d*d)/length(d)))
    df_tmp=data.frame(den$x,den$y)
    colnames(df_tmp)=c(paste0("c_",n,"_x"),paste0("c_",n,"_y"))
    df=cbind(df,df_tmp)
    #lines(density(sqrt(colSums(d*d)/length(d))),col=colors[col])
    #fig <- fig %>% add_trace(density(sqrt(colSums(d*d)/length(d))),y = ~y,x=~x, name = 'trace 0',mode = 'lines')
    col=col+1
  }
  legend("topright",legend = names(data_var),title=v,lty=1,col=colors[c(1:col)])
}

fig=plot_ly(df,type="scatter")
for (n in names(data_var)){
  fig <- fig %>% add_trace(x = df[[paste0("c_",n,"_x")]], y= df[[paste0("c_",n,"_y")]], name = n, mode = 'lines')
}

fig2=plot_ly(df,type="scatter")
for (n in names(data_var)){
  fig2 <- fig2 %>% add_trace(x = df[[paste0("c_",n,"_x")]], y= df[[paste0("c_",n,"_y")]], name = n, mode = 'lines')
}
fig <- subplot(fig, fig2)
fig

names(data_var)



lines(density(sqrt(colSums(d*d)/length(d))))
lines(den$x,den$y)
fig <- plot_ly(x=den$x,y=den$y, name = 'trace 0',mode = 'lines')
fig <- fig %>% add_trace(x=den$x,y=den$y, name = 'trace 1',mode = 'scatter',col="red")
fig


par(mar=c(3,3.5,2,0.5),mfrow=c(2,2))
for (v in variables){
  data_var=summarize_by_var_age(age,v,1:650,init)
  plot(age_borehole[,c(2,1)],type="l",col=1,xlab="",ylab="",main="Age")
  col=1
  for (n in names(data_var)){
    d=data_var[[n]]
    polygon(c(d[,which(colSums(d)==min(colSums(d)))],
              rev(d[,which(colSums(d)==max(colSums(d)))])),
            c(5*(1:650),5*rev(1:650)),col=colors[col],border=NA)
    col=col+1
  }
  legend("topright",legend = names(data_var),title=v,fill=colors[c(1:col)])
}
library(plotly)
plot_ly(x=c(1:650),y=data_var[[n]][,1],type="scatter")

plot(c(1:650),data_var[[n]][,1])

data_var[[n]][,1]


xlength(data_var[[n]][,1])

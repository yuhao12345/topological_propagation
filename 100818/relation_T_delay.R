library(R.matlab)
library('dplyr')
library(hexbin)

df<-readMat('E:/dwell/0/27_0.mat')
df1<-as.data.frame(df)

df2<-df1 %>%
  filter(temp.1>0)

writeMat(con='C:/Users/ykang/Dropbox/TI7/T_delay.mat', data=df2)
# bin<-hexbin(df2$temp.1,df2$temp.2, xbins=250)
# 
# plot(bin,xlab="Time delay", ylab="Transmittance")

plot(df2$temp.1,df2$temp.2,xlab="Time delay", ylab="Transmittance",
col=rgb(0,0,100,25,maxColorValue=255), pch=16)

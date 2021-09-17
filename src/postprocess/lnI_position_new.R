# lnI   has spin up and down

library(R.matlab)
library('dplyr')
library(doParallel)
library(foreach)

# length<-30
# for nnn, spin up/down half is 0, for rashba, all are nonzero
path_test <- "E:/dwell2/48/" # path to data files

cores=detectCores()
cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl)


df<-readMat(paste0(path_test,0,'.mat'))
dfu<-as.data.frame(df$u)
df1<-dfu %>%     #get position
  # filter(abs(V1)<length) %>%
  filter(V2>10) %>%
  group_by(V1)  %>%
  summarize_all(sum)

## integral over corss section first, then log, then do average for all channels in all samples

dos  <-  foreach (m = 0:1999,.combine=cbind,.packages=c('dplyr','R.matlab')) %dopar% {
  df<-readMat(paste0(path_test,m,'.mat'))
  dfu<-as.data.frame(df$u)
  dfd<-as.data.frame(df$d)
  # head(df1)
  dfu[,3:ncol(dfu)]<-dfu[,3:ncol(dfu)]^2
  df1<-dfu %>%     #get position
    # filter(abs(V1)<length) %>%
    # filter(V2>10) %>%
    #arrange(V1)  %>%
    group_by(V1)  %>%
    summarize_all(sum)
  ln_u<-rowMeans (log(df1[,3:12]),na.rm=TRUE)   #log  3:ncol(df1)  3:12
  dfd[,3:ncol(dfd)]<-dfd[,3:ncol(dfd)]^2
  df2<-dfd %>%     #get position
    # filter(abs(V1)<length) %>%
    filter(V2>10) %>%
    group_by(V1)  %>%
    summarize_all(sum)
  ln_d<-rowMeans (log(df2[,13:22]),na.rm=TRUE)  #3:ncol(df1)   13:22
  u<-as.data.frame(df1$V1)
  u$lnu<-ln_u
  u$lnd<-ln_d
  u[,2:3]
  # write.csv(u,paste0(path_test,m,'_0.csv')) 
  
}


#stop cluster
stopCluster(cl)
f<-cbind(df1$V1,dos) 
write.csv(f, 'E:/dwell2/3_8.csv')    

# logI: number or number_1 , I (y>0): number_0  
# I(y=disorder region): number_3   I(y>width-7):number_4   I(y>width-10):number_5 
# I(y>width-8):number_6   logI(y>width-8):number_7 
# logI(y>width-5): number_8
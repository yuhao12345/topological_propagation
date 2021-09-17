# lnI   no spin

library(R.matlab)
library('dplyr')
library(doParallel)
library(foreach)

# length<-30
# for nnn, spin up/down half is 0, for rashba, all are nonzero
path_test <- "E:/dwell2/26/" # path to data files

# path_test <- 'E:/narrowlead/11/'
cores=detectCores()
cl <- makeCluster(8) #not to overload your computer
registerDoParallel(cl)


df<-readMat(paste0(path_test,0,'.mat'))
dfu<-as.data.frame(df$u)
df_tmp<-dfu %>%     #get position
  # filter(abs(V1)<length) %>%
  filter(V2>0) %>%
  group_by(V1)  %>%
  summarize_all(sum)

## integral over corss section first, then log, then do average for all channels in all samples


dos  <-  foreach (m = 0:1999,.combine=cbind,.packages=c('dplyr','R.matlab')) %dopar% {
  df<-readMat(paste0(path_test,m,'.mat'))
  dfu<-as.data.frame(df$u)

  dfu[,3:ncol(dfu)]<-(dfu[,3:ncol(dfu)]^2)
  df1<-dfu %>%     #get position
    filter(V2>0) %>%
    group_by(V1)  %>%
    summarize_all(sum)
  ln_u<-rowMeans ((df1[,3:ncol(dfu)]),na.rm=TRUE)   #log  3:ncol(df1)
  u<-as.data.frame(df_tmp$V1)
  u$lnu<-ln_u

  u[,2]
  
}

#stop cluster
stopCluster(cl)
f<-cbind(df_tmp$V1,dos)
write.csv(f, 'E:/dwell2/26_0.csv')

# write.csv(f, 'E:/narrowlead/17_0.csv') 
# logI: number  I: number_0 
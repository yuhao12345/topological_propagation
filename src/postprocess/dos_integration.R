# dos

library(R.matlab)
library('dplyr')
library(doParallel)
library(foreach)

length<-30

path_test <- "E:/dwell/11/" # path to data files

cores=detectCores()
cl <- makeCluster(5) #not to overload your computer
registerDoParallel(cl)



## integral over corss section first, then log, then do average for all channels in all samples


dos  <-  foreach (m = 0:4999,.combine=cbind,.packages=c('dplyr','R.matlab')) %dopar% {
  df<-readMat(paste0(path_test,m,'.mat'))
  df<-as.data.frame(df)
  # head(df1)
  
  df1<-df %>%     #get position
    filter(abs(ans.1)<length) %>%
    filter(ans.2>0) #%>%
    # group_by(ans.1)  %>% 
    # summarize(sum=sum(ans.3))
  # df1$sum
  sum(df1$ans.3)   #dos within -length~length
}

#stop cluster
stopCluster(cl)

write.csv(dos, 'E:/dwell/11.csv')
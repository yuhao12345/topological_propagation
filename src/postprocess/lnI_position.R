#run C:\Users\user\Dropbox\kwant\100818\QVH_QSH_QVH_new.py  first

library(R.matlab)
library('dplyr')
library(doParallel)
library(foreach)

length<-180

path_test <- "E:/yuhao/65/" # path to data files

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# lnt_m <- data.frame()

## log at each point, then average
# for (m in 0:14){
#   df<-readMat(paste0('C:/Users/ykang/Documents/yuhao_share/30/',m,'.mat'))
#   
#   df1<-as.data.frame(df)
#   # head(df1)
#   
#   lnt<-df1 %>% 
#     filter(abs(ans.1)<length) %>%
#     # filter(ans.2<=15) %>%
#     arrange(ans.1)  %>%
#     group_by(ans.1) %>% 
#     summarize(mean=mean(ans.3))
#   if (m==0){
#     lnt_m<-lnt
#   }
#   else{
#     lnt_m <- merge(lnt_m, lnt,by="ans.1") 
#   }
#   
# }

## integral over corss section first, then log, then do average for all channels in all samples
df0<-readMat(paste0(path_test,'0.mat'))
df<-as.data.frame(df0)
# head(df1)

pos<-df %>%     #get position
  filter(abs(ans.1)<length) %>%
  # filter(ans.2<=15) %>%
  arrange(ans.1)  %>%
  group_by(ans.1)  %>% 
  summarize(sum=sum(ans.3))

lnt_m  <-  foreach (m = 0:15999,.combine=cbind,.packages=c('dplyr','R.matlab')) %dopar% {
  df<-readMat(paste0(path_test,m,'.mat'))
  
  df1<-as.data.frame(df)
  # head(df1)
  
  inte<-df1 %>%     #integration over cross scetion
    filter(abs(ans.1)<length) %>%
    # filter(ans.2<=15) %>%
    arrange(ans.1)  %>%
    group_by(ans.1)  %>% 
    summarize(sum=sum(abs(ans.3)^2))
  
  inte$sum
}

#stop cluster
stopCluster(cl)
lnt_m<-cbind(pos$ans.1,lnt_m)
# for (m in 0:2){
#   df<-readMat(paste0('C:/Users/ykang/Documents/yuhao_share/30/',m,'.mat'))
#   
#   df1<-as.data.frame(df)
#   # head(df1)
#   
#   inte<-df1 %>%     #integration over cross scetion
#     filter(abs(ans.1)<length) %>%
#     # filter(ans.2<=15) %>%
#     arrange(ans.1)  %>%
#     group_by(ans.1)  %>% 
#     summarize(sum=sum(abs(ans.3)^2))
#   if (m==0){
#     lnt_m<-inte
#   }
#   else{
#     lnt_m <- merge(lnt_m, inte,by="ans.1") 
#   }
#   
# }
 write.csv(lnt_m, 'E:/yuhao/lnt65.csv')


# plot(lnt$ans.1,lnt$mean)


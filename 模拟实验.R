get_test_res <- function(par){
  #par=c(miu=1,A=1,B=-1.2,ia=0.8,ib=-0.5,sd=1)
  n_samples=par[1]
  n_snps=par[2]
simulated_data=get_simdata(n_samples=par[1],n_snps=par[2],par=par[-c(1:2)])
res_single<-sapply(1:n_snps,function(x)single_snp_res(simulated_data[,x],pheno = simulated_data$Y))

res_double<-sapply(1:2,function(x)double_snp_res(x,simulated_data))
thre <- 0.05
  
mode1_power <- length(which(!is.na(match(1:2,which(res_single[2,] < thre))
)))
mode1_fpr <- length(which(!is.na(match(3:5,which(res_single[2,] < thre))
)))

mode2_power <- length(which(!is.na(match(2,which(res_double[,1] < thre))
)))+length(which(!is.na(match(3,which(res_double[,2] < thre))
)))
mode2_fpr <- length(which(!is.na(match(c(1:5)[-2],which(res_double[,1] < thre))
)))+length(which(!is.na(match(c(1:5)[-3],which(res_double[,2] < thre))
)))
return(c(mode1_power,mode1_fpr,mode2_power,mode2_fpr))


}


library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=100,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.05)),cl=cl)
stopCluster(cl)

n100_sd0.05_res <-  rowSums(res)/(c(2,3,2,6)*100)




library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=100,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.1)),cl=cl)
stopCluster(cl)

n100_sd0.1_res <-  rowSums(res)/(c(2,3,2,6)*100)



library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=100,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.5)),cl=cl)
stopCluster(cl)

n100_sd0.5_res <-  rowSums(res)/(c(2,3,2,6)*100)





library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=300,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.05)),cl=cl)
stopCluster(cl)

n300_sd0.05_res <-  rowSums(res)/(c(2,3,2,6)*100)




library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=300,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.1)),cl=cl)
stopCluster(cl)

n300_sd0.1_res <-  rowSums(res)/(c(2,3,2,6)*100)



library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=300,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.5)),cl=cl)
stopCluster(cl)

n300_sd0.5_res <-  rowSums(res)/(c(2,3,2,6)*100)





library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=600,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.05)),cl=cl)
stopCluster(cl)

n600_sd0.05_res <-  rowSums(res)/(c(2,3,2,6)*100)




library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=600,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.1)),cl=cl)
stopCluster(cl)

n600_sd0.1_res <-  rowSums(res)/(c(2,3,2,6)*100)



library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=600,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.5)),cl=cl)
stopCluster(cl)

n600_sd0.5_res <-  rowSums(res)/(c(2,3,2,6)*100)







library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=1000,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.05)),cl=cl)
stopCluster(cl)

n1000_sd0.05_res <-  rowSums(res)/(c(2,3,2,6)*100)




library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=1000,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.1)),cl=cl)
stopCluster(cl)

n1000_sd0.1_res <-  rowSums(res)/(c(2,3,2,6)*100)



library(parallel)
cl.cores <- 12
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("lihood",'lihood0',"get_simdata",'double_interaction_res','single_snp_res','get_test_res',"or_res","double_snp_res"))
clusterEvalQ(cl,{library(mvtnorm)})
clusterEvalQ(cl,{library(tidyverse)})

res <- pbsapply(1:100,function(c)get_test_res(par=c(n_samples=1000,n_snps=5,miu=0.3,A=0.8,B=-0.9,ia=0.6,ib=-0.5, sd=0.5)),cl=cl)
stopCluster(cl)

n1000_sd0.5_res <-  rowSums(res)/(c(2,3,2,6)*100)




dat_sd0.05 <-   rbind(n100_sd0.05_res,n300_sd0.05_res,n600_sd0.05_res,n1000_sd0.05_res)

dat_sd0.1 <- rbind(n100_sd0.1_res,n300_sd0.1_res,n600_sd0.1_res,n1000_sd0.1_res)
dat_sd0.5 <- rbind(n100_sd0.5_res,n300_sd0.5_res,n600_sd0.5_res,n1000_sd0.5_res)



dat_power <- data.frame(sample_size=c(100,300,600,1000),variance_0.05=dat_sd0.05[,3],variance_0.1=dat_sd0.1[,3],variance_0.5=dat_sd0.5[,3])



dat_fpr <- data.frame(sample_size=c(100,300,600,1000),variance_0.05=dat_sd0.05[,4],variance_0.1=dat_sd0.1[,4],variance_0.5=dat_sd0.5[,4])


library(ggplot2)
library(tidyr)
library(dplyr)

df_long <- dat_power %>%
  pivot_longer(cols = starts_with("variance"), 
               names_to = "variance", 
               values_to = "power")

pdf("power_plot.pdf",width = 7,height = 6)
ggplot(df_long, aes(x = sample_size, y = power, color = variance, group = variance)) +
  geom_line() + 
  geom_point(size=2) +
  labs(x = "Sample Size", y = "Power", color = "Error Variance") +xlab(label = NULL)+ylab(label = NULL)+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    legend.position = "top"
  )+ theme_bw() +
  scale_x_continuous(breaks = c(0,100,300,600,1000))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()
  )+theme(axis.text.x = element_text( color="black", size=20),axis.text.y =element_text( color="black", size=20),
          panel.border = element_rect(color = "black", size = 1.5, fill = NA),
          axis.ticks.length.y = unit(0.2,"cm"),axis.ticks.length.x = unit(0.2,"cm") ,axis.title.x = element_text(size = 20),  # 修改横坐标标题的字体大小
          axis.title.y = element_text(size = 20))

dev.off()






df_long <- dat_fpr %>%
  pivot_longer(cols = starts_with("variance"), 
               names_to = "variance", 
               values_to = "FPR")

pdf("fpr_plot.pdf",width = 7,height = 6)
ggplot(df_long, aes(x = sample_size, y = FPR, color = variance, group = variance)) +
  geom_line() + 
  geom_point(size=2) +
  labs(x = "Sample Size", y = "FPR", color = "Error Variance") +xlab(label = NULL)+ylab(label = NULL)+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    legend.position = "top"
  )+ theme_bw() +
  scale_y_continuous(position = "right") +
  scale_x_continuous(breaks = c(0,100,300,600,1000))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()
  )+theme(axis.text.x = element_text( color="black", size=20),axis.text.y =element_text( color="black", size=20),
          panel.border = element_rect(color = "black", size = 1.5, fill = NA),
          axis.ticks.length.y = unit(0.2,"cm"),axis.ticks.length.x = unit(0.2,"cm") ,axis.title.x = element_text(size = 20),  # 修改横坐标标题的字体大小
          axis.title.y = element_text(size = 20))

dev.off()



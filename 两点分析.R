library(parallel)
library(pbapply)
library(mvtnorm)
 double_interaction_res<- function(pheno_dat,snpdat){



    a00 <- which(snpdat[,1]==0&snpdat[,2]==0)######genotype
    a01 <- which(snpdat[,1]==0&snpdat[,2]==1)
    a10 <- which(snpdat[,1]==1&snpdat[,2]==0)
    a11 <- which(snpdat[,1]==1&snpdat[,2]==1)

    ##########################################aa
    c_aa <- rep(NA,length(pheno_dat))
    c_aa[c(a00,a11)] <- 1/2
    c_aa[c(a01,a10)] <- -1/2
    
    
    
    dat <- cbind(pheno_dat,c_aa)
    o_r <-  or_res(dat)
    
    parres_h1<-  optim(par=c(0.5,0.5),lihood,y=dat[,1],c=dat[,2],
                            
                            method = "BFGS",control = list(maxit=20000)
      )
     
    parres_h0<-  optim(par=c(0.5),lihood0,y=dat[,1],
                           
                            method = "BFGS",control = list(maxit=20000)
      )
      
      l1 <- parres_h1$value
      l0 <- parres_h0$value
      lr <- 2*(l0-l1)
    if(lr<0){
      lr <- 0
    }
    pvalue <- pchisq(lr,df=1,lower.tail = F)

    # res<- c(
    #   lr=lr,pvalue=pvalue,or_min=o_r[1],or=o_r[2],or_max=o_r[3])
    # 
    res<- c(
      pvalue)
    

   
    

  return(res) 
}

 double_snp_res<- function(wpos,simulated_data){
   #wpos=2
   pheno_dat=simulated_data$Y
   interaction_snp = simulated_data[,-c(wpos,ncol(simulated_data))]
   
   res<- sapply(1:(ncol(simulated_data)-2),function(c){
     
    snpdat = cbind(simulated_data[,wpos],interaction_snp[,c])
     res = double_interaction_res(pheno_dat,snpdat)

   })
   
  return( res)

   
 }
   
 double_snp_res(1,simulated_data)
 double_snp_res(2,simulated_data)
 
   
   

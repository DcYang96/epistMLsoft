lihood <- function(y,A,c){

  fenzi <- exp(A[1]+A[2]*c
               );
  fenmu <- 1 + fenzi;
  py1 <- fenzi/fenmu;
  py0 <- 1 - py1; 
  
  id1 <- which(y == 1);
  id0 <- which(y == 0);
  logL <- sum(log(py1[id1])) + sum(log(py0[id0]))
  
  return(-logL)}
lihood0 <- function(y,A){
  # A=0.5
  # y=dat[,1]
  fenzi <- exp(A[1]
               );
  fenmu <- 1 + fenzi;
  py1 <- rep(fenzi/fenmu,length(y))
  py0 <- rep(1 - py1,length(y))
  
  id1 <- which(y == 1);
  id0 <- which(y == 0);
  logL <- sum(log(py1[id1])) + sum(log(py0[id0]))
  
  return(-logL)}
or_res <- function(or_matrix){
  
  odd_ratio <- c(a=length(which(or_matrix[,1]==1&or_matrix[,2]==1/2)),b=length(which(or_matrix[,1]==1&or_matrix[,2]==-1/2)),
                 c=length(which(or_matrix[,1]==0&or_matrix[,2]==1/2)),d=length(which(or_matrix[,1]==0&or_matrix[,2]==-1/2))
  )
  odd_ratio[which(odd_ratio==0)]=1
  sqrt(sum(1/odd_ratio))
  o_r <- (odd_ratio[1]*odd_ratio[4])/(odd_ratio[2]*odd_ratio[3])
  or_min <- exp(log(o_r)-1.96* sqrt(sum(1/odd_ratio)))
  or_max <- exp(log(o_r)+1.96* sqrt(sum(1/odd_ratio)))
  return(c(or_min,o_r,or_max)) 
}
single_snp_res <- function(snp,pheno){
  # snp=simulated_data[,1]
  # pheno = simulated_data$Y
  a0 <- which(snp==0)######genotype
  a1 <- which(snp==1)
 

  
  dat <- cbind(pheno,snp)
  dat[which(dat[,2]==1),2] = 1/2
  dat[which(dat[,2]==0),2] = -1/2
  o_r <-  or_res(dat)
  # dat[,9] <-  dat[,9]/2
 parres_h1<-  optim(par=c(0.5,0.5),lihood,y=dat[,1],c=dat[,2],method = "BFGS",control = list(maxit=20000)
  )
  # optim(par=rep(0,9),lihood,y=dat[,1],c=rep(0,nrow(dat)),x1=dat[,2]
  # ,x2=dat[,3],x3=dat[,4],x4=dat[,5],x5=dat[,6],x6=dat[,7]
  #,x7=dat[,8],method = "BFGS",control = list(maxit=20000)
  #)
 parres_h0<-  optim(par=c(0.5),lihood0,y=dat[,1],method = "BFGS",control = list(maxit=20000)
  )
  
 l1 <- parres_h1$value
 l0 <- parres_h0$value
 lr <- 2*(l0-l1)
 pvalue <- pchisq(lr,df=1,lower.tail = F)

  res<- c(
    lr=lr,pvalue=pvalue,or_min=o_r[1],or=o_r[2],or_max=o_r[3])
  
  return(res)
}
res_single<-sapply(1:5,function(c)single_snp_res(simulated_data[,c],pheno = simulated_data$Y))

res_single

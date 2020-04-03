####spherical GMM(sigma=1)
d=3
n=100
k=2
mean_t=cbind(c(0,0,0),c(2,3,4))
prob_t=c(.3,.7)
#common cov=1
sigma_t1=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3)
sigma_t2=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3)
label=seq(k)
index=sample(label,size = n,replace = TRUE,prob = prob_t)
require(mvtnorm)
x=rmvnorm(n,mean = mean_t[,1],sigma = sigma_t1)
x[which(index==2),]=rmvnorm(sum(index==2),mean=mean_t[,2],sigma = sigma_t2)

####initial values
mean_0=cbind(colMeans(x[which(index==1),]),colMeans(x[which(index==2),]))
prob_0=rep(1/k,k)
sigma_1=cov(x[which(index==1),])
sigma_2=cov(x[which(index==2),])
sigma_0=list(sigma_1,sigma_2)
####EM Algorithm
GMM_EM=function(x,mean_0,sigma_0,prob_0){
  it=0
  error=0.01
  converge=FALSE
  n=nrow(x)
  k=2
  d=ncol(x)
  mean=mean_0
  sigma=sigma_0
  prob=prob_0
  log_l=0
  weight=matrix(0,ncol=k,nrow=n)
  showit=TRUE
  
  if (showit){
    cat(paste("iterations of EM","\n"))
  }
  
  while (!converge ){
    for (j in 1:k){
      weight[,j]=prob[j]*dmvnorm(x,mean = mean[,j],sigma = sigma[[j]])
      
    }
    weight=weight/rowSums(weight)
    prob=colSums(weight)/n
    for (j in 1:k){
      mean[,j]=(t(x)%*%(weight[,j]))/sum(weight[,j])
      var=matrix(0,ncol=d,nrow=d)
      for (i in 1:n){
        var=var+weight[i,j]*(x[i,]-mean[,j])%*%t(x[i,]-mean[,j])
      }
      sigma[[j]]=var/sum(weight[,j])
    }
    it=it+1
    log_old=log_l
    for (i in 1:n){
      a=NULL
      for (j in 1:k){
        a[j]=prob[j]*dmvnorm(x[i,],mean[,j],sigma[[j]],log = FALSE)
      }
      
      log_l=log_l+log(sum(a))
    }
    
    if (abs(log_l-log_old)<=error){
      converge=TRUE
    }
  }
  cluster=rep(1,n)
  for (i in i:n){
    cluster[i]=which.max(weight[i,])
  }
  return(list(mean=mean,sigma=sigma,prob=prob,cluster=cluster,it=it))
}
###speed of EM
start_EM=Sys.time()
params_EM=GMM_EM(x,mean_0,sigma_0,prob_0)
end_EM=Sys.time()
time_delta_EM=end-start
##accuracy of EM
accuracy=mean(params_EM$cluster==index)
m1=abs(params_EM$sigma[[1]]-diag(rep(1,3)))
m2=abs(params_EM$sigma[[2]]-diag(rep(1,3)))
error_max_EM=max(m1,m2,abs(params_EM$mean-mean_t),abs(params_EM$prob-prob_t))


#####tenser moments estimate###############
GMM_TM=function(x){
  k=2
  d=dim(x)[2]
  n=dim(x)[1]
  x_half1=x[1:(n/2),]
  x_half2=x[(n/2+1):n,]
  mu_hat=colMeans(x_half1)
  M2_hat=(t(x_half1)%*%(x_half1))/nrow(x_half1)
  mat=M2_hat-mu_hat%*%t(mu_hat)
  var_hat=eigen(mat)$values[2]
  svd=svd(mat)
  m2_hat=svd$u[,1:2]%*%diag(svd$d[1:2])%*%t(svd$v[,1:2])
  U_hat=svd(m2_hat)$u[,1:2]
  library(MASS)
  m=t(U_hat)%*%m2_hat%*%U_hat
  m_sqrt=eigen(m)$vectors%*%diag(sqrt(eigen(m)$values))%*%t(eigen(m)$vectors)
  B_hat=U_hat%*%m_sqrt
  W_hat=U_hat%*%ginv(m_sqrt)
  mu_hat_2=colMeans(x_half2)
  c=t(W_hat)%*%mu_hat_2
  ##Define the third moment
  M3=function(x,y,z,dim_array){
    array=array(0,dim=c(dim_array[1],dim_array[2],dim_array[3]))
    for (i in 1:dim_array[1]){
      for (j in 1:dim_array[2]){
        for (k in 1:dim_array[3]){
          array[i,j,k]=x[i]*y[j]*z[k]
        }
      }
    }
    return(array)
  }
  
  ##
  M_3_hat=array(0,dim=c(k,k,k))
  for (i in 1:dim(x_half2)[1]){
    x_w=t(W_hat)%*%x_half2[i,]
    dim_m=rep(dim(x_w)[1],3)
    M_3_hat=M_3_hat+M3(x_w,x_w,x_w,dim_m)
  }
  M_3_hat=M_3_hat/(n/2)
  ##
  e1=c(1,0,0)
  e2=c(0,1,0)
  e3=c(0,0,1)
  w1=t(W_hat)%*%e1
  w2=t(W_hat)%*%e2
  w3=t(W_hat)%*%e3
  w=list(w1,w2,w3)
  s=array(0,dim=c(k,k,k))
  for (i in 1:d){
    dim_s=rep(2,3)
    s=s+M3(c,w[[i]],w[[i]],dim_s)+M3(w[[i]],c,w[[i]],dim_s)+M3(w[[i]],w[[i]],c,dim_s)
  }
  m_3_hat=M_3_hat-var_hat*s
  ###Iterative process
  delta=0.05
  key=NULL
  params=list()
  for (t in 1:ceiling(log(1/delta,base = 2))) {
    theta_1=runif(1,min=-1,max=1)
    theta_2=runif(1,-sqrt(1-theta_1^2),sqrt(1-theta_1^2))
    theta=c(theta_1,theta_2)
    
    M_3_theta=array(0,dim=c(2,2,1))
    for (i in 1:dim(x_half2)[1]){
      x_w=t(W_hat)%*%x_half2[i,]
      x_t=t(W_hat%*%theta)%*%x_half2[i,]
      M_3_theta=M_3_theta+M3(x_w,x_w,x_t,c(2,2,1))
    }
    M_3_theta=M_3_theta/(n/2)
    
    w1_theta=t(W_hat%*%theta)%*%e1
    w2_theta=t(W_hat%*%theta)%*%e2
    w3_theta=t(W_hat%*%theta)%*%e3
    w_theta=list(w1_theta,w2_theta,w3_theta)
    s_theta=array(0,dim=c(2,2,1))
    for (i in 1:d){
      s_theta=s_theta+M3(c,w[[i]],w_theta[[i]],c(2,2,1))+M3(w[[i]],c,w_theta[[i]],c(2,2,1))+M3(w[[i]],w[[i]],t(W_hat%*%theta)%*%mu_hat_2,c(2,2,1))
    }
    m_3_theta=M_3_theta-var_hat*s_theta
    m_3_theta=as.matrix(m_3_theta[,,1])
    e_values=eigen(m_3_theta)$values
    e_vectors=as.matrix(eigen(m_3_theta)$vectors)
    r=length(e_values)
    set=NULL
    for (i in 1:r){
      for (j in 1:r){
        set[(r-1)*i+j]=abs(e_values[i]-e_values[j])
      }
    }
    
    key[t]=which.min(c(set,abs(e_values)))
    mu_hat_theta=matrix(0,nrow=d,ncol=k)
    for (i in 1:k){
      mu_hat_theta[,i]=(e_values[i]*(B_hat%*%e_vectors[,i]))/(t(theta)%*%e_vectors[,i])
    }
    w_hat_theta=ginv(mu_hat_theta)%*%mu_hat
    params[[t]]=list(mu_hat_theta,w_hat_theta)
  }
  index_best=which.max(key)
  params_best=params[[index_best]]
  return(params_best)
}
##speed of TM
start_TM=Sys.time()
params_TM=GMM_TM(x)
end_TM=Sys.time()
time_delta_TM=end_TM-start_TM
##accuracy of TM
params_true=list(mean_t,prob_t)
cov_true=1
error_max_TM=max(abs(params_TM-params_true),abs(cov_true-var_hat))
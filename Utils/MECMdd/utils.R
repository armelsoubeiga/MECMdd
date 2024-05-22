library(R.utils)

#' Creation of a matrix of focal sets
makeF<- function(c,type=c('simple','full','pairs'),pairs=NULL,Omega=TRUE){
  if(type=='full'){
    ii<-1:2^c
    N<-length(ii)
    F<-matrix(0,N,c)
    CC<-intToBin(0:(N-1))
    for(i in 1:N) F[i,]<-as.numeric(substring(CC[i],1:c,1:c))
    F<-F[,c:1]
  }else{
    F<-rbind(rep(0,c),diag(c))
    if(type=='pairs'){
      if(is.null(pairs)){
        for(i in 1:(c-1)){
          for(j in (i+1):c){
            f<-rep(0,c)
            f[c(i,j)]<-1
            F<-rbind(F,f)
          }
        }
      } else{
        n<-nrow(pairs)
        for(i in 1:n){
          f<-rep(0,c)
          f[pairs[i,]]<-1
          F<-rbind(F,f)
        }
      }
    }
    if(Omega & !((type=="pairs")&(c==2)) & !((type=="simple")&(c==1))){
      F<-rbind(F,rep(1,c))
    }
  }
  row.names(F)<-NULL
  return(F)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
extractMass <- function(mass,F,g=NULL,S=NULL,method,crit=NULL,Kmat=NULL,trace=NULL,
                        D=NULL,W=NULL,J=NULL,param=NULL){
  
  n<-nrow(mass)
  c<-ncol(F)
  if(any(F[1,]==1)){
    F<-rbind(rep(0,c),F) 
    mass<-cbind(rep(0,n),mass)
  }
  f<-nrow(F)
  card<-rowSums(F)
  
  conf<-mass[,1]            
  C<- 1/(1-conf)
  mass.n<- C*mass[,2:f]   
  pl<- mass%*% F          
  pl.n<- C*pl             
  p <-pl/rowSums(pl)      
  bel<- mass[,card==1]   
  bel.n<-C*bel            
  y.pl<-max.col(pl)       
  y.bel<-max.col(bel)     
  Y<-F[max.col(mass),]    
  # non dominated elements
  Ynd<-matrix(0,n,c)
  for (i in 1:n){
    ii<-which(pl[i,]>= bel[i,y.bel[i]])
    Ynd[i,ii]<-1
  }
  
  P<-F/card
  P[1,]<- rep(0,c)
  betp<- mass %*% P   
  betp.n<- C* betp    
  
  lower.approx<-vector(mode='list',length=c)
  upper.approx<-vector(mode='list',length=c)
  lower.approx.nd<-vector(mode='list',length=c)
  upper.approx.nd<-vector(mode='list',length=c)
  nclus<-rowSums(Y)
  outlier<-which(nclus==0) # outliers
  nclus.nd<-rowSums(Ynd)
  for(i in 1:c){
    upper.approx[[i]]<- which(Y[,i]==1)                 
    lower.approx[[i]]<- which((Y[,i]==1) & (nclus==1))  
    upper.approx.nd[[i]]<- which(Ynd[,i]==1)                 
    lower.approx.nd[[i]]<- which((Ynd[,i]==1) & (nclus.nd==1)) 
  }
  # Nonspecificity
  card<-c(c,card[2:f])
  Card<-matrix(card,n,f,byrow=TRUE)
  N <- sum(log(Card)*mass)/log(c)/n
  clus<-list(conf=conf,F=F,mass=mass,mass.n=mass.n,pl=pl,pl.n=pl.n,bel=bel,bel.n=bel.n,y.pl=y.pl,
             y.bel=y.bel,Y=Y,betp=betp,betp.n=betp.n,p=p,upper.approx=upper.approx,
             lower.approx=lower.approx,Ynd=Ynd,upper.approx.nd=upper.approx.nd,
             lower.approx.nd=lower.approx.nd,N=N,outlier=outlier,g=g,S=S,
             crit=crit,Kmat=Kmat,trace=trace,D=D,method=method,W=W,J=J,param=param)
  class(clus)<-"credpart"
  return(clus)
}



#%%%%%%%%%%%
lambdaInit_local <- function(weight, p, c){
  lambda <- matrix(0, nrow = c, ncol = p)
  for(j in 1:c){
    if(weight=='sum'){
      lambda[j,] <- rep((1/p),p)
    }else if(weight=='prod'){
      lambda[j,] <- rep(1,p)
    }
  }
  rownames(lambda) <- paste0("k_",1:c)
  colnames(lambda) <- paste0("P_",1:p)
  return(lambda)
}


#%%%%%%%%%%%
lambdaInit_local <- function(weight, p, c){
  lambda <- matrix(0, nrow = c, ncol = p)
  for(j in 1:c){
    if(weight=='sum'){
      lambda[j,] <- rep((1/p),p)
    }else if(weight=='prod'){
      lambda[j,] <- rep(1,p)
    }
  }
  rownames(lambda) <- paste0("k_",1:c)
  colnames(lambda) <- paste0("P_",1:p)
  return(lambda)
}


lambdaInit_global <- function(weight,p){
  lambda <- numeric(p)
  
  if(weight=='sum'){
    lambda <- rep((1/p),p)
  }else if(weight=='prod'){
    lambda <- rep(1,p)
  }
  return(lambda)
}



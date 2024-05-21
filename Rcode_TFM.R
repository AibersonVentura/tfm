install.packages("VGAM")
library(VGAM)
## Variables ##

T=5 #minimum length for the estimation period that banks are allowed to use
B=10*10^6 #Number of simulations
N = 50*10^3 # Portfolio Size
w=0.24 # correlation. The maximum value provided by the Regulation
PDLR<- c(0.3, 1, 5, 10)/100 #Long-run Probability of Default
alpha<- c(99,99.5,99.9)/100 # Confidence level Alpha
DRt<-c() # Default rates at each time period
vec_q<-matrix(0,nrow=B, ncol=length(alpha))
DR_var<-DR_mean <- matrix(0, nrow=B, ncol=length(PDLR))
WCDR<- matrix(0, nrow=length(alpha), ncol=length(PDLR))
rownames(WCDR)<- paste("Alpha", 100*alpha, sep=" ")
colnames(WCDR)<- paste("PDLR", 100*PDLR, sep=" ")
WCDR2<-WCDR

## WCDR with unknown PDLR ##

# Importance of Sampling 
set.seed(2023)
for(i in 1:length(PDLR)){
  for (g in 1:B) {
    zt <- rnorm(5) # Sample from standard normal distribution
    fz <- pnorm((qnorm(PDLR[i]) - sqrt(w) * zt) / sqrt(1 - w))  
    s=qnorm(PDLR[i])
    Theta= max(0,((s)*(1-fz))/(fz*(N-s)))
    if(Theta!=0) fz<- (fz*exp(Theta))/(1+fz*(exp(Theta)-1))  # Conditional default probability
    DRz <- vapply(fz, function(x) rbinom(1, N, x), as.integer(1L))/N  # Default rate at t conditioned on Z
    DRt <- mean(DRz) # Estimated PDLR
    
    DR_mean[g,i] <- DRt  # Mean for the UpperBound
    DR_var[g,i]<- var(DRz)# Variance for the UpperBound
    
    for(k in 1:length(alpha)){ 
      vec_q[g,k]<-pnorm((qnorm(DRt)-sqrt(w)*qnorm(1-alpha[k]))/sqrt(1-w))
    }
    
  }
  
  WCDR[,i]<- colMeans(vec_q)
}

WCDR

# Antitéticas

DR_var2<-DR_mean2 <- matrix(0, nrow=B, ncol=length(PDLR))
empieza <- Sys.time()
set.seed(2023)
for(i in 1:length(PDLR)){
  
  for (g in 1:B) {
    zt <- rnorm(5) # Sample from standard normal distribution
    fz1 <- pnorm((qnorm(PDLR[i]) - sqrt(w) * zt) / sqrt(1 - w))
    fz2<- pnorm((qnorm(PDLR[i]) - sqrt(w) * (-zt)) / sqrt(1 - w))
    fz<-  (fz1+fz2)/2  # Conditional default probability
    DRz <- vapply(fz, function(x) rbinom(1, N, x), as.integer(1L))/N  # Default rate at time t conditioned on Z
    DRt <- mean(DRz) # Estimated PDLR
    DR_mean2[g,i] <- DRt  #UB
    DR_var2[g,i]<- var(DRz)#UB
    for(k in 1:length(alpha)){ 
      vec_q[g,k]<-pnorm((qnorm(DRt)-sqrt(w)*qnorm(1-alpha[k]))/sqrt(1-w))
    }
    
  }
  
  WCDR2[,i]<- colMeans(vec_q)
}

WCDR2

## WCDR with known PDLR ##
alpha<-c(95,99,99.5,99.9)/100
q<- matrix(0, ncol=length(PDLR), nrow=length(alpha))

for(i in 1:length(alpha)){
  
  for(j in 1:length(PDLR)){
    
    q[i,j]<- pnorm((qnorm(PDLR[j])-sqrt(w)*qnorm(1-alpha[i]))/sqrt(1-w))
  }
}

colnames(q)<- c("PDLR 0.3%","PDLR 1%", "PDLR 5%", "PDLR 10%" )
rownames(q)<- c("Alpha 95%","Alpha 99%", "Alpha 99.5%", "Alpha 99.9%")
round(q*100,2)

####

## Upper Bound
#DR_mean es una matriz donde cada columna son los diferentes DRt conseguidos con Montecarlo para una Probabilidad de default PDLR determinada
#load("C:/Users/aiber/OneDrive/Desktop/TFM_articles/data.Rdata")
alpha<-c(95,99,99.5,99.9)/100
Upbounds <- list()
beta<-seq(0.51,0.57, by=0.002)
q_b<- matrix(0, ncol=length(PDLR), nrow=length(alpha))

colnames(q_b)<-paste("PDLR", PDLR, sep=" ")
rownames(q_b)<-paste("Alpha", alpha, sep=" ")

q_es_up<- matrix(0, ncol=length(beta), nrow=length(alpha))
colnames(q_es_up)<-paste("Beta", beta, sep=" ")
rownames(q_es_up)<-paste("Alpha", alpha, sep=" ")

for (k in 1:length(PDLR)) {
  
  s<- qnorm(PDLR[k]) 
  
  varianza<-(pbinorm(s, s, mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = w) - pnorm(s)^2)
  
  for(i in 1:(length(alpha))){
    
    for(j in 1:length(beta)){
      
      UPB <- DR_mean[,k] + qnorm(beta[j])*sqrt(varianza/T)
      qs<- pnorm((qnorm(UPB)-sqrt(w)*qnorm(1-alpha[i]))/sqrt(1-w))
      q_es_up[i,j]<- mean(qs)
      
      if(q_es_up[i,j]>=q[i,k] & q_b[i,k]==0) q_b[i,k]<-beta[j]
    }
    
  } 
  
  Upbounds[[k]]<- q_es_up  
}
names(Upbounds)<- paste("PDLR", PDLR, sep=" ")
Upbounds ## Quantiles for each Beta.

q_b # Smallest beta >= real value.


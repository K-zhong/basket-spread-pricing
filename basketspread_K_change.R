rm(list = ls())

library(BB)
library(nleqslv)
library(openxlsx)

### case1----
n = 1
m = 2

s1 = c(95)
s2 = c(90)
s3 = c(105)
corr = matrix(c(1, 0.9, 0.8,
                0.9, 1, 0.9,
                0.8, 0.9, 1), nrow = 3, byrow = TRUE)

r = 0.03
T = 1
t = 0

w = matrix(c(1,
             0.8,
             0.2), nrow = 3, byrow = TRUE) 

S= matrix(c(s1, s2, s3), nrow = 3, byrow = TRUE) 

sig = 0.2245*matrix(1,nrow = 3, ncol =1) 

cov_mat = corr * (sig %*% t(sig))
rw =  matrix(c(1,
               -0.8,
               -0.2), nrow = 3, byrow = TRUE) 

Ib = colSums(rw*S)


MT = seq(0.75, 1.2, by = 0.025)

KK = matrix(MT) %*%t(matrix(Ib))

del = matrix(c(-0.9, -0.7578, 0.10, 0.7578, 0.90))

### share----

kappa = sqrt(2*T/pi)


cli = function(sigi,deltai)
{ 
  A = log(2* pnorm(sigi*deltai*sqrt(T), mean = 0, sd = 1) )
}

I_owen = function(om, osig , oa, oa1, oa2)
{ 
  om = as.vector(om)
  osig = as.vector(osig)
  oa = as.vector(oa)
  oa1 = as.vector(oa1)
  oa2 = as.vector(oa2)
  inerf = function(u){ exp( oa*u - (u - om)^2/(2*osig^2) ) * pnorm(oa1+oa2*u, mean = 0, sd = 1) }
  result = (integrate( inerf, lower = 0, upper = Inf)$value ) / (osig*sqrt(2*pi))
}

### Theorem 5----

cov_mat1 = cov_mat[1:n,1:n] 
fcov2 = function(m){
  if(m==0) {A = c()}
  else{ A = cov_mat[(n+1):(n+m),(n+1):(n+m)] }
  return(A)
}

cov_mat2 = fcov2(m)

upsilon1 = matrix(c(rep(1, times = n)), nrow = 1, ncol = n)
upsilon2  = matrix(c(rep(-1, times = m)), nrow = 1, ncol = m)
upsilonk = c(upsilon1,upsilon2)

phiB = function(x ){
  A = 2*exp(0.5*x^2*T)* pnorm(sqrt(T)*x, mean = 0, sd = 1)
}

L = function(m,a,delta,bk,CK){
  if (m==0) {A= phiB(0) /( a *exp(0.5*sum( ( bk[c(1:n)]* sqrt(1-delta[1:n]^2) )%*%t(bk[c(1:n)]*sqrt(1-delta[1:n]^2)) * cov_mat1 *T ) ) *phiB(sum(bk[c(1:n)]* sig[1:n]*delta[1:n] ) )   ) }
  else{ A=(exp(0.5*sum( ( bk[c((n+1):(n+m))]* sqrt(1-delta[(n+1):(n+m)]^2) )%*%t(bk[c((n+1):(n+m))]*sqrt(1-delta[(n+1):(n+m)]^2)) * cov_mat2 *T ) ) *phiB(sum(bk[c((n+1):(n+m))]* sig[(n+1):(n+m)]*delta[(n+1):(n+m)] ) )  )/( a *exp(0.5*sum( ( bk[c(1:n)]* sqrt(1-delta[1:n]^2) )%*%t(bk[c(1:n)]*sqrt(1-delta[1:n]^2)) * cov_mat1 *T ) ) *phiB(sum(bk[c(1:n)]* sig[1:n]*delta[1:n] ) )   ) 
  }
  return(A)
}

sig_LB = function(m,delta,bk,CK){
  if ( m==0){A =  sqrt(sum( ( upsilonk[c(1:n)]*bk[c(1:n)]* sqrt(1-delta[1:n]^2) )%*%t( upsilonk[c(1:n)]*bk[c(1:n)]*sqrt(1-delta[1:n]^2)) * cov_mat1 *T ) ) }  
  else{
    A = sqrt( sum( ( upsilonk[c(1:(n+m))]*bk[c(1:(n+m))]* sqrt(1-delta^2) )%*%t( upsilonk[c(1:(n+m))]*bk[c(1:(n+m))]*sqrt(1-delta^2)) * cov_mat *T ) )
  }
  return(A)
}

fetaup = function(m,delta,upsilonk,bk,CK){
  if(m==0){A = sum(upsilonk[1:n]*bk[1:n]*sig[1:n]*delta[1:n])/sig_LB(m,delta,bk,CK)}
  else{A = sum(upsilonk[1:(n+m)]*bk[1:(n+m)]*sig[1:(n+m)]*delta[1:(n+m)])/sig_LB(m,delta,bk,CK) }
  return(A)
}

fb2 = function(m){
  if(m==0){A = c() }
  else{A = w[(n+1):(n+m)]* S[(n+1):(n+m)]*exp(r*T)/( w[(n+1):(n+m)]%*% S[(n+1):(n+m)] *exp(r*T)  + CK)}
  return(A)
}

fa = function(m,CK){
  if(m==0){A =   CK/ ( w[1:n]%*% S[1:n] *exp(r*T)  ) }
  else{A =  ( w[(n+1):(n+m)]%*% S[(n+1):(n+m)] *exp(r*T)  + CK)/ ( w[1:n]%*% S[1:n] *exp(r*T)  ) }
  return(A)
}

Dih = function(i,u,delta)
{
  A = S[i] * exp( (r- 0.5*sig[i]^2*delta[i]^2)*T - cli(sig[i],delta[i]) + sig[i]*delta[i]*u  )
}

d0h = function(u,a,delta,bk,CK)
{
  A =  ( log(L(m,a,delta,bk,CK))+  sum(upsilonk[c(1:(n+m))]*bk[c(1:(n+m)) ]*sig*delta )*u )/ sig_LB(m,delta,bk,CK)
}

dih = function(i,u,a,delta,bk,CK)
{
  A =  ( log(L(m,a,delta,bk,CK))+  sum(upsilonk[c(1:(n+m))]*bk[c(1:(n+m)) ]*sig*delta )*u + sum(upsilonk[c(1:(n+m))]*bk[c(1:(n+m))]*corr[,i]*sig*sqrt(1-delta^2)*T ) *sig[i]*sqrt(1-delta[i]^2) )/ sig_LB(m,delta,bk,CK)
}

result_arr_bas1_thrp_LB = array(0, dim = c(length(MT), length(del), 2 ) )
dimnames(result_arr_bas1_thrp_LB) <- list( MT, del, c("LB","CLB") )

for (l in 1 : length(del) ) {
  delta = matrix(c( rep(del[l], times = (n+m) ) ) )
  bc_l = matrix(cli(sig,delta))
  for (j in 1 : length(KK) ) {
    
    CK = KK[j]
    b1 = w[1:n]* S[1:n]*exp(r*T)/( w[1:n]%*% S[1:n] *exp(r*T) )
    b2 = fb2(m)
    bk = c(b1,b2)
    etaup = fetaup(m,delta,upsilonk,bk,CK)
    a = fa(m,CK) 
    Vi_LB = c()
    for (i in 1 : length(S) ) { 
      Vi_LB[i] = Dih(i ,0 , delta )*2* I_owen(0, sqrt(T) , sig[i]*delta[i], dih( i, 0,a,delta,bk,CK), etaup) 
    }
    V0_LB = 2* I_owen(0, sqrt(T) , 0, d0h(0,a,delta,bk,CK), etaup) 
    LB = exp(-r*T) * (t(w*upsilonk)%*%Vi_LB - CK*V0_LB )
    VF = t(w*upsilonk)%*%S - CK*exp(-r*T)
    result_arr_bas1_thrp_LB[j,l,"LB"] = LB
    CLB = max(LB,0) + max( VF- LB,0)
    result_arr_bas1_thrp_LB[j,l,"CLB"] = CLB
  }
}

### Theorem 5 with Remark 2----
result_arr_bas1_thrp_MLB = array(0, dim = c(length(MT), length(del), 1 ) )
dimnames(result_arr_bas1_thrp_MLB) <- list( MT, del, c("MLB") )

cov_mat1 = cov_mat[1:n,1:n] 
fcov2 = function(m){
  if(m==0) {A = c()}
  else{ A = cov_mat[(n+1):(n+m),(n+1):(n+m)] }
  return(A)
}

cov_mat2 = fcov2(m)

upsilon1 = matrix(c(rep(1, times = n)), nrow = 1, ncol = n)
upsilon2  = matrix(c(rep(-1, times = m)), nrow = 1, ncol = m)
upsilonk = c(upsilon1,upsilon2)

phiB = function(x ){
  A = 2*exp(0.5*x^2*T)* pnorm(sqrt(T)*x, mean = 0, sd = 1)
}

L = function(m,a,delta,bk,CK){
  if (m==0) {A= phiB(0) /( a *exp(0.5*sum( ( bk[c(1:n)]* sqrt(1-delta[1:n]^2) )%*%t(bk[c(1:n)]*sqrt(1-delta[1:n]^2)) * cov_mat1 *T ) ) *phiB(sum(bk[c(1:n)]* sig[1:n]*delta[1:n] ) )   ) }
  else{ A=(exp(0.5*sum( ( bk[c((n+1):(n+m))]* sqrt(1-delta[(n+1):(n+m)]^2) )%*%t(bk[c((n+1):(n+m))]*sqrt(1-delta[(n+1):(n+m)]^2)) * cov_mat2 *T ) ) *phiB(sum(bk[c((n+1):(n+m))]* sig[(n+1):(n+m)]*delta[(n+1):(n+m)] ) )  )/( a *exp(0.5*sum( ( bk[c(1:n)]* sqrt(1-delta[1:n]^2) )%*%t(bk[c(1:n)]*sqrt(1-delta[1:n]^2)) * cov_mat1 *T ) ) *phiB(sum(bk[c(1:n)]* sig[1:n]*delta[1:n] ) )   ) 
  }
  return(A)
}

sig_LB = function(m,delta,bk,CK){
  if ( m==0){A =  sqrt(sum( ( upsilonk[c(1:n)]*bk[c(1:n)]* sqrt(1-delta[1:n]^2) )%*%t( upsilonk[c(1:n)]*bk[c(1:n)]*sqrt(1-delta[1:n]^2)) * cov_mat1 *T ) ) }  
  else{
    A = sqrt( sum( ( upsilonk[c(1:(n+m))]*bk[c(1:(n+m))]* sqrt(1-delta^2) )%*%t( upsilonk[c(1:(n+m))]*bk[c(1:(n+m))]*sqrt(1-delta^2)) * cov_mat *T ) )
  }
  return(A)
}

fetaup = function(m,delta,upsilonk,bk,CK){
  if(m==0){A = sum(upsilonk[1:n]*bk[1:n]*sig[1:n]*delta[1:n])/sig_LB(m,delta,bk,CK)}
  else{A = sum(upsilonk[1:(n+m)]*bk[1:(n+m)]*sig[1:(n+m)]*delta[1:(n+m)])/sig_LB(m,delta,bk,CK) }
  return(A)
}

fb2 = function(m){
  if(m==0){A = c() }
  else{A = w[(n+1):(n+m)]* S[(n+1):(n+m)]*exp(r*T)/( w[(n+1):(n+m)]%*% S[(n+1):(n+m)] *exp(r*T)  + CK)}
  return(A)
}

fa = function(m,CK){
  if(m==0){A =   CK/ ( w[1:n]%*% S[1:n] *exp(r*T)  ) }
  else{A =  ( w[(n+1):(n+m)]%*% S[(n+1):(n+m)] *exp(r*T)  + CK)/ ( w[1:n]%*% S[1:n] *exp(r*T)  ) }
  return(A)
}

Dih = function(i,u,delta)
{
  A = S[i] * exp( (r- 0.5*sig[i]^2*delta[i]^2)*T - cli(sig[i],delta[i]) + sig[i]*delta[i]*u  )
}


md0h = function(u,a,x,delta,bk,CK)
{
  A =  ( x +  sum(upsilonk[c(1:(n+m))]*bk[c(1:(n+m)) ]*sig*delta )*u )/ sig_LB(m,delta,bk,CK)
}

mdih = function(i,u,a,x,delta,bk,CK)
{
  A =  ( x +  sum(upsilonk[c(1:(n+m))]*bk[c(1:(n+m)) ]*sig*delta )*u + sum(upsilonk[c(1:(n+m))]*bk[c(1:(n+m))]*corr[,i]*sig*sqrt(1-delta^2)*T ) *sig[i]*sqrt(1-delta[i]^2) )/ sig_LB(m,delta,bk,CK)
}

for (l in 1 : length(del) ) {
  delta = matrix(c( rep(del[l], times = (n+m) ) ) )
  bc_l = matrix(cli(sig,delta))
  for (j in 1 : length(KK) ) {

    CK = KK[j]
    b1 = w[1:n]* S[1:n]*exp(r*T)/( w[1:n]%*% S[1:n] *exp(r*T) )
    b2 = fb2(m)
    bk = c(b1,b2)
    etaup = fetaup(m,delta,upsilonk,bk,CK)
    a = fa(m,CK) 
    Vi_LB = c()
    VF = t(w*upsilonk)%*%S - CK*exp(-r*T)
    fun = function(x) {
      for (i in 1 : length(S) ) { 
        Vi_LB[i] = Dih(i ,0 , delta )*2* I_owen(0, sqrt(T) , sig[i]*delta[i], mdih( i, 0,a,x,delta,bk,CK), etaup) 
      }
      V0_LB = 2* I_owen(0, sqrt(T) , 0, md0h(0,a,x,delta,bk,CK), etaup) 
      LB = exp(-r*T) * (t(w*upsilonk)%*%Vi_LB - CK*V0_LB )
      result = max(LB,0) + max( VF- LB,0)
      return(-result)
    }
    
    result = BBoptim(par = c(log(L(m,a,delta,bk,CK))), fn = fun, lower = c(-Inf), upper = c(Inf))
    MLB = -result$value
    result_arr_bas1_thrp_MLB[j,l,"MLB"] = MLB
  }
}



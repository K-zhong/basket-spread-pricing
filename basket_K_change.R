rm(list = ls())

library(BB)
library(nleqslv)
library(openxlsx)

### Kchange_1----
n= 4
m = 0
corr = 0.75*matrix(1,nrow = n, ncol =n) 
diag(corr) = 1
r = 0.03
T = 1
t = 0

w = 0.25*matrix(1,nrow = n, ncol =1) 
S = 100*matrix(1,nrow = n, ncol =1) 
sig = 0.2245*matrix(1,nrow = n, ncol =1)  
cov_mat = corr * (sig %*% t(sig))
rw = 0.25*matrix(1,nrow = n, ncol =1) 

Ib = colSums(rw*S)
del = matrix(c(0.7578))

MT = seq(0.5, 1.5, by = 0.025)

KK = matrix(MT) %*%t(matrix(Ib)) 


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
    #K = CK[j]
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



### share function----
sig_lambda = function(delta)
{ 
  A = sqrt( sum( (w* sqrt(1-delta^2) )%*%t(w*sqrt(1-delta^2)) * cov_mat *T ))
}

rho_ilam = function(i,delta)
{
  A =  t( w* sqrt(1-delta^2) * sig )%*%corr[,i]/sig_lambda(delta)*sqrt(T)
}

sig_ilam = function(i,delta)
{
  if ( numToInts(rho_ilam(i,delta))[1] == 1 ) {
    A = 0 } else {A = sqrt( ( 1- rho_ilam(i,delta)^2 )*T )}
}

Di = function(i,u,delta)
{
  A = S[i] * exp( (r- 0.5*sig[i]^2)*T - cli(sig[i],delta[i]) + sig[i]*delta[i]*u + 0.5*sig[i]^2*(1-delta[i]^2)*sig_ilam(i,delta)^2 + 0.5*sig[i]^2*rho_ilam(i,delta)^2*(1-delta[i]^2)*T  )
}

d0 = function(u,v,delta)
{
  A = - ( v - t(w) %*% (sig*delta)* (u- kappa) )/ sig_lambda(delta)
}

di = function(i,u,v,delta)
{
  A = - ( v - rho_ilam(i,delta)*sig[i]*sig_lambda(delta)*sqrt(1-delta[i]^2)*sqrt(T) - t(w) %*% (sig*delta)*(u- kappa)  )/ sig_lambda(delta)
}


### Theorem 1----
result_arr_bas1_thrp_I1 = array(0, dim = c(length(MT), length(del), 1) )
dimnames(result_arr_bas1_thrp_I1) <- list( MT, del, c("I1") )

for (l in 1 : length(del) ) {
 
  delta = matrix(c( rep(del[l], times = (n+m) ) ) )
  
  bc_l = matrix(cli(sig,delta))
  
  eta = t(w) %*% (sig*delta) / sig_lambda(delta)
  
  for (j in 1 : length(KK) ) {

    CK = KK[j]
    d_lambda = log(CK) - t(w) %*% ( log(S) + (r-0.5*sig^2)*T - bc_l + kappa*sig*delta) #1个数值 
   
    Vi_1 = c()
    for (i in 1 : length(S) ) { 
      Vi_1[i] = Di(i ,0 , delta )*2* I_owen(0, sqrt(T) , sig[i]*delta[i], di( i, 0,d_lambda,delta), eta) 
    }
    
    V0_1 = 2* I_owen(0, sqrt(T) , 0 , d0( 0,d_lambda,delta), eta) 
   
    I1 = exp(-r*T) * (t(w)%*%Vi_1 - CK*V0_1 )
    
    result_arr_bas1_thrp_I1[j,l,"I1"] = I1
  }
}

### Theorem 3----
result_arr_bas1_thrp_mIC = array(0, dim = c(length(MT), length(del), 1) )
dimnames(result_arr_bas1_thrp_mIC) <- list( MT, del, c("mIC") )

phiB = function(x ){
  A = 2*exp(0.5*x^2*T)* pnorm(sqrt(T)*x, mean = 0, sd = 1)
}

fic = function(i,u,kappa,delta,eta)
{
  A =  S[i] * exp( (r- 0.5*sig[i]^2)*T - cli(sig[i],delta[i]) + 0.5*sig[i]^2*(1-delta[i]^2)*sig_ilam(i,delta)^2 + sig[i]*rho_ilam(i,delta)*sqrt(1-delta[i]^2)*sqrt(T) *( u/sig_lambda(delta) + kappa* eta )  ) *phiB(  sig[i]*delta[i] - eta*sig[i]*sqrt(1-delta[i]^2)*rho_ilam(i,delta)*sqrt(T) ) 
}


for (l in 1 : length(del) ) {

  delta = matrix(c( rep(del[l], times = (n+m) ) ) )
  
  bc_l = matrix(cli(sig,delta))
  eta = t(w) %*% (sig*delta) / sig_lambda(delta)
  
  for (j in 1 : length(KK) ) {
 
    CK = KK[j]
    ffic = function(x) {
      t(w)%*%matrix( c(fic(c(1),x,kappa,delta,eta), fic(c(2),x,kappa,delta,eta) , fic(c(3),x,kappa,delta,eta), fic(c(4),x,kappa,delta,eta)) , nrow = (n+m), byrow = TRUE) - CK
    }
    xstartp = c(0)
    resultp = nleqslv(xstartp, ffic)
    lambp = resultp$x
    
    Vi_2mC = c()
    for (i in 1 : length(S) ) {
      Vi_2mC[i] = Di(i ,0 , delta )*2* I_owen(0, sqrt(T) , sig[i]*delta[i], di( i, 0,lambp,delta), eta)
    }
    
    V0_2mC = 2* I_owen(0, sqrt(T) , 0 , d0( 0,lambp,delta), eta)
    
    mIC = exp(-r*T) * (t(w)%*%Vi_2mC - CK*V0_2mC )
    
    result_arr_bas1_thrp_mIC[j,l,"mIC"] = mIC
  }
}


### Theorem 4----
result_arr_bas1_thrp_IA = array(0, dim = c(length(MT), length(del), 1) )
dimnames(result_arr_bas1_thrp_IA) <- list( MT, del, c("IA") )

fpsi = function(i,u,kappa,delta,eta)
{
  A =  S[i] * sig[i]*sqrt(1-delta[i]^2)*rho_ilam(i,delta)*sqrt(T)*exp( (r- 0.5*sig[i]^2)*T - cli(sig[i],delta[i]) + 0.5*sig[i]^2*(1-delta[i]^2)*sig_ilam(i,delta)^2 + sig[i]*rho_ilam(i,delta)*sqrt(1-delta[i]^2)*sqrt(T) *( u/sig_lambda(delta) + kappa* eta ) + 0.5*( sig[i]*delta[i] - eta*sig[i]*sqrt(1-delta[i]^2)*rho_ilam(i,delta)*sqrt(T) )^2*T ) * 2* pnorm( ( sig[i]*delta[i] - eta*sig[i]*sqrt(1-delta[i]^2)*rho_ilam(i,delta)*sqrt(T) )*sqrt(T), mean = 0, sd = 1)
}

flambi = function(i,u,kappa,delta,eta)
{
  A =  S[i] * exp( (r- 0.5*sig[i]^2)*T - cli(sig[i],delta[i]) + 0.5*sig[i]^2*(1-delta[i]^2)*sig_ilam(i,delta)^2 + sig[i]*rho_ilam(i,delta)*sqrt(1-delta[i]^2)*sqrt(T) *( u/sig_lambda(delta) + kappa* eta ) + 0.5*( sig[i]*delta[i] - eta*sig[i]*sqrt(1-delta[i]^2)*rho_ilam(i,delta)*sqrt(T) )^2*T ) * 2* pnorm( ( sig[i]*delta[i] - eta*sig[i]*sqrt(1-delta[i]^2)*rho_ilam(i,delta)*sqrt(T) )*sqrt(T), mean = 0, sd = 1)
}

for (l in 1 : length(del) ) {
  
  delta = matrix(c( rep(del[l], times = (n+m) ) ) )
  
  bc_l = matrix(cli(sig,delta))
  
  eta = t(w) %*% (sig*delta) / sig_lambda(delta)
  ff = function(x) {
    t(w)%*%matrix( c(fpsi(c(1),x,kappa,delta,eta), fpsi(c(2),x,kappa,delta,eta), fpsi(c(3),x,kappa,delta,eta), fpsi(c(4),x,kappa,delta,eta)  ) , nrow = 4, byrow = TRUE)
  }
  xstart = c(0)
  result = nleqslv(xstart, ff)
  
  fv =  t(w)%*% matrix( c(flambi(c(1),result$x,kappa,delta,eta), flambi(c(2),result$x,kappa,delta,eta) ,flambi(c(3),result$x,kappa,delta,eta) ,flambi(c(4),result$x,kappa,delta,eta)  ) , nrow = 4, byrow = TRUE)
  
  Vi_2A = c()
  for (j in 1 : length(KK) ) {
    CK = KK[j]
    if ( fv >= CK ) {
      for (i in 1 : length(S) ) { Vi_2A[i] =  S[i] }
      V0_2C = 1
    } else {
     
      ffv = function(x) {
        x = c(x[1], x[2])
        f1 = t(w)%*%matrix( c(flambi(c(1),x[1],kappa,delta,eta), flambi(c(2),x[1],kappa,delta,eta),flambi(c(3),x[1],kappa,delta,eta),flambi(c(4),x[1],kappa,delta,eta)   ) , nrow = 4, byrow = TRUE) - CK
       
        f2 = t(w)%*%matrix( c(flambi(c(1),x[2],kappa,delta,eta), flambi(c(2),x[2],kappa,delta,eta),flambi(c(3),x[2],kappa,delta,eta),flambi(c(4),x[2],kappa,delta,eta)  ) , nrow = 4, byrow = TRUE) - CK
     
        c(f1,f2)
      }
      xstartv = c(-5,5)
      resultv = nleqslv(xstartv, ffv)
      d_lamb1 = resultv$x[1]
      d_lamb2 = resultv$x[2]
      
      if (abs(d_lamb1 - d_lamb2)  <= 0.0000001 ){
       
        for (i in 1 : length(S) ) {
          Vi_2A[i] = Di(i ,0 , delta )*2* I_owen(0, sqrt(T) , sig[i]*delta[i], di( i, 0,d_lamb2,delta), eta)
        }
        V0_2C = 2*  I_owen(0, sqrt(T) , 0, d0( 0,d_lamb2,delta), eta)
      } else {
        for (i in 1 : length(S) ) {
          Vi_2A[i] = Di(i ,0 , delta )*2* ( I_owen(0, sqrt(T) , sig[i]*delta[i], -di( i, 0,d_lamb1,delta), -eta) +  I_owen(0, sqrt(T) , sig[i]*delta[i], di( i, 0,d_lamb2,delta), eta)  )
        }
        V0_2C = 2* ( I_owen(0, sqrt(T) , 0, -d0(  0,d_lamb1,delta), -eta) +  I_owen(0, sqrt(T) , 0, d0( 0,d_lamb2,delta), eta)  )
      }
      
    }
    IA =  exp(-r*T) * (t(w)%*%Vi_2A - CK*V0_2C )
  
    result_arr_bas1_thrp_IA[j,l,"IA"] = IA
  }
  
}


### Theorem 2----
result_arr_bas1_thrp_IM = array(0, dim = c(length(MT), length(del), 1) )
dimnames(result_arr_bas1_thrp_IM) <- list( MT, del, c("IM") )

ht = function(x)
{
  A =  exp(0.5*x^2*T)* 2 * pnorm(x*sqrt(T), mean = 0, sd = 1)
}

ft = function(x,delta)
{
  A =  exp( x *t(w)%*% ( log(S) + (r- 0.5*sig^2)*T - cli(sig,delta) ) + 0.5*x^2* sum( (w* sqrt(1-delta^2) )%*%t(w*sqrt(1-delta^2)) * cov_mat *T ) )*ht( t(w) %*% (sig*delta) *x)
}

DM = function(u,c,m,delta)
{
  A = exp( c* t(w)%*% ( log(S) + (r- 0.5*sig^2)*T - cli(sig,delta) +  kappa* sig*delta ) + m + 0.5*c^2*sig_lambda(delta)^2 - c*(kappa - u) * t(w) %*% (sig*delta) )
}

dM = function(u,v,c,delta)
{
  A = ( v - c*sig_lambda(delta)^2 - t(w) %*% (sig*delta)*(u- kappa)  )/ sig_lambda(delta)
}

d0M = function(u,v,delta)
{
  A = ( v  - t(w) %*% (sig*delta)*(u- kappa)  )/ sig_lambda(delta)
}

M1 = t(w)%*%S*exp(r*T)
M2 =  function(delta)
{
  A = 0
  for (i in 1:length(S)){
    for (j in 1:length(S)){
      A = A + w[i]*w[j]*S[i]*S[j]*exp( T*(2*r  - 0.5*sig[i]^2*delta[i]^2 -  0.5*sig[j]^2*delta[j]^2 + corr[i,j]*sig[i]*sig[j]*sqrt(1-delta[i]^2)*sqrt(1-delta[j]^2)  -cli(sig[i],delta[i]) -cli(sig[j],delta[j])   )  ) *ht(sig[i]*delta[i]+sig[j]*delta[j])
    }
  }
  return(A)
}
M3 = function(delta)
{
  A = 0
  for (i in 1:length(S)){
    for (j in 1:length(S)){
      for (k in 1: length(S)) {
        A = A + w[i]*w[j]*w[k]*S[i]*S[j]*S[k]*exp( T*(3*r  - 0.5*sig[i]^2*delta[i]^2 - 0.5*sig[j]^2*delta[j]^2 -  0.5*sig[k]^2*delta[k]^2 + corr[i,j]*sig[i]*sig[j]*sqrt(1-delta[i]^2)*sqrt(1-delta[j]^2) + corr[i,k]*sig[i]*sig[k]*sqrt(1-delta[i]^2)*sqrt(1-delta[k]^2) + corr[k,j]*sig[k]*sig[j]*sqrt(1-delta[k]^2)*sqrt(1-delta[j]^2) -cli(sig[i],delta[i]) -cli(sig[j],delta[j]) -cli(sig[k],delta[k])  )  ) *ht(sig[i]*delta[i]+sig[j]*delta[j]+sig[k]*delta[k])
      }
    }
  }
  return (A)
}


for (l in 1 : length(del) ) {
  delta = matrix(c( rep(del[l], times = (n+m) ) ) )
  
  MM2 = M2(delta)
  MM3 = M3(delta)
  eta_ms = (MM3-3*M1*MM2+ 2*M1^3)/(MM2 - M1^2)^(1.5)
  fc = function(x) {
    ft(3*x,delta) + 2*ft(x,delta)^3 - 3*ft(x,delta)*ft(2*x,delta) - eta_ms*(ft(2*x,delta) - ft(x,delta)^2)^(1.5)
  }
  
  cstart = c(2)
  cresult = nleqslv(cstart, fc)
  c = cresult$x
  mm = 0.5*log((MM2 - M1^2)/(ft(2*c,delta) - ft(c,delta)^2 )  )
  tau = M1 - exp(mm)*ft(c,delta)
  
  bc_l = matrix(cli(sig,delta))
  
  etam = c*t(w) %*% (sig*delta)
  eta = t(w) %*% (sig*delta) / sig_lambda(delta)
  
  for (j in 1 : length(KK) ) {
    #K = CK[j]
    CK = KK[j]
    d_lambda = log(CK) - t(w) %*% ( log(S) + (r-0.5*sig^2)*T - bc_l + kappa*sig*delta) 
    d_hlambda = log(CK-tau)/c- mm/c - t(w) %*% ( log(S) + (r-0.5*sig^2)*T - bc_l + kappa*sig*delta) 
    
    if ( c >= 0 ) {
      if (d_lambda >= d_hlambda)
      {V_2M = DM(0,c,mm,delta)* 2*( I_owen( 0, sqrt(T), etam, dM(0,d_lambda,c,delta), -eta) -  I_owen( 0, sqrt(T), etam, dM(0,d_hlambda,c,delta), -eta) )
      V0_2M = 2*( I_owen( 0, sqrt(T), 0, d0M(0,d_lambda,delta), -eta) -  I_owen( 0, sqrt(T), 0, d0M(0,d_hlambda,delta), -eta) ) }
      else{ V_2M =0
      V0_2M =0 }
    } else {
      if (d_lambda >= d_hlambda){ V_2M = DM(0,c,mm,delta)* 2*I_owen( 0, sqrt(T), etam, dM(0,d_hlambda,c,delta), -eta)
      V0_2M = 2*I_owen( 0, sqrt(T), 0, d0M(0,d_hlambda,delta), -eta) }
      else{ V_2M = DM(0,c,mm,delta)* 2* I_owen( 0, sqrt(T), etam, dM(0,d_lambda,c,delta), -eta)
      V0_2M = 2* I_owen( 0, sqrt(T), 0, d0M(0,d_lambda,delta), -eta) }
    }
    
    IM = exp(-r*T) * (V_2M + (tau- CK)*V0_2M ) + result_arr_bas1_thrp_I1[j,l,]
   
    result_arr_bas1_thrp_IM[j,l,"IM"] = IM
  }
}




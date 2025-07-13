clear all; clc;
% Case 3
n = 10^8; d = 3;
n2 = 10^5;
filename1 = 'basketspread_case3_nonMC.xlsx';
data = readtable(filename1);
LB_data = table2array(data(:,2));

sig = 0.2245*ones(1,d); 
S = [100, 90, 95];

corr = [1, 0.9, 0.8;
        0.9, 1, 0.9;
        0.8, 0.9, 1]; 

cov_mat = diag(sig)*corr*diag(sig);

r = 0.03;
del = 0.7578;
T = 1;
rw = [0.6, 0.4, -1];
w = [0.6, 0.4, 1]; 
K = [0.5:0.1:1.5]'*(rw*S'); 
np = sum(rw>0);
nn = sum(rw<0);
sigi = sig(1);

kappa = sqrt(2*T/pi);

bas1_se = zeros(size(K,1),size(del,2)); %每一行是一个篮子在不同K下的se
bas1_price = zeros(size(K,1),size(del,2));
bas1_seLB = zeros(size(K,1),size(del,2)); %每一行是一个篮子在不同K下的se
bas1_priceLB = zeros(size(K,1),size(del,2));

phibt = @(u) 2*exp(0.5*u^2*T)* normcdf(sqrt(T)*u, 0, 1);

rate = zeros(2,size(K,1));

for l = 1:size(del,2)

    delta = [repmat(del(l),1,d)]';%[del(l), del(l), del(l)]';
    deltai = delta(1);
    li = log(2* normcdf(sigi*deltai*sqrt(T),  0,  1) ) ;

    D_delta = sqrt( eye(d) - diag(delta)^2 ); 
 
    ncorr_mat = [corr, zeros(1,d)';zeros(1,d),1];
    nmu = zeros(1,d+1);
    nrv1 = randn(n, d+1);  
    L = chol(ncorr_mat, 'lower');  
    nrv1 = nmu + nrv1 * L'; 
    
    rv1 = nrv1(:,1:d);
    rv2 = nrv1(:,d+1);

    nrv12 = randn(n2, d+1);
    nrv12 = nmu + nrv12 * L'; 
    rv12 = nrv12(:,1:d);
    rv22 = nrv12(:,d+1);

    X = ( D_delta*rv1' + delta*abs(rv2') )' ;
    X2 = ( D_delta*rv12' + delta*abs(rv22') )' ;
     
    cl = log( 2* normcdf(sig'.*delta*sqrt(T)) ); 
    coef = log(S') + (r-0.5*sig'.*sig')*T - cl;
    ST =  exp(coef + repmat(sig' ,1,n).*X'*sqrt(T)); 
    ST2 =  exp(coef + repmat(sig' ,1,n2).*X2'*sqrt(T));

    term5 = w(1:np).*S(1:np)*exp( (r - 0.5*sigi^2)*T - li );
    term52 = w(end- nn+1 :end).*S(end- nn+1 :end)*exp( (r - 0.5*sigi^2)*T - li );


    for i = 1:length(K)
        LB = LB_data(i);
        payoff2 = exp(-r*T).*max(0, rw*ST2- K(i)); 
       
        b1 = w(1:np).*S(1:np)*exp(r*T)/(w(1:np)*(S(1:np)')*exp(r*T) );
        b2 = w(end- nn+1 :end).*S(end- nn+1 :end)*exp(r*T)/(w(end- nn+1 :end)*(S(end- nn+1 :end)')*exp(r*T)+K(i) );
        vecb = [b1, b2]; 
        L1_term5 = prod(term5.^vecb(1:np),2)*exp(0.5*T*(1-deltai^2)*( sum(sum(vecb(1:np)'*vecb(1:np).*cov_mat(1:np,1:np))))) * phibt( vecb(1:np)*sig(1:np)'*deltai );
        L1_term52 = prod(term52.^vecb(end- nn+1 :end),2)*exp(0.5*T*(1-deltai^2)*( sum(sum(vecb(end- nn+1 :end)'*vecb(end- nn+1 :end).*cov_mat(end- nn+1 :end,end- nn+1 :end))))) * phibt( vecb(end- nn+1 :end)*sig(end- nn+1 :end)'*deltai );
        LA_term5 =  prod(ST(1:np,:).^(vecb(1:np)'), 1)*prod(w(1:np).^(vecb(1:np)), 2); 
        LA_term6 =  prod(ST(end- nn+1 :end,:).^(vecb(end- nn+1 :end)'), 1)*prod(w(end- nn+1 :end).^(vecb(end- nn+1 :end)), 2); 
        LA = LA_term5./LA_term6; 

        LA2_term5 =  prod(ST2(1:np,:).^(vecb(1:np)'), 1)*prod(w(1:np).^(vecb(1:np)), 2); 
        LA2_term6 =  prod(ST2(end- nn+1 :end,:).^(vecb(end- nn+1 :end)'), 1)*prod(w(end- nn+1 :end).^(vecb(end- nn+1 :end)), 2); 
        LA2 = LA2_term5./LA2_term6; 

        aa = ( w(end- nn+1 :end)*(S(end- nn+1 :end)'*exp(r*T)) + K(i)) /( w(1:np)*(S(1:np)')*exp(r*T));
        L1 = aa * L1_term5/L1_term52;
        contro_LB = (rw*ST2- K(i)).*( LA2 >= (L1) ) *exp(-r*T);
        covLB = cov(contro_LB, payoff2);
        rate(1,i) = -covLB(1,2)/var(contro_LB);
        
        payoff_LB = exp(-r*T).*( max(0, rw*ST- K(i)) + rate(1,i)*(rw*ST- K(i)).*( LA >= (L1)  )) - rate(1,i)*LB; 
        [cp_LB, cv_LB, ci_LB] = normfit(payoff_LB);
        bas1_priceLB(i,l) = cp_LB;
        bas1_seLB(i,l) = (ci_LB(2)-ci_LB(1))/cp_LB; 

        payoff = exp(-r*T).*max(0, rw*ST- K(i)); 
        [cp, cv, ci] = normfit(payoff);
        bas1_price(i,l) = cp;
        bas1_se(i,l) = (ci(2)-ci(1))/cp;

    end
end

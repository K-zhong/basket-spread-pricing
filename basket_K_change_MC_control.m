clear all; clc;

n = 10^8; d = 4;
n2 = 10^5;
sig = 0.2245*ones(1,d); 
S = 100*ones(1,d); 
corr = 0.75 *ones(d);
corr(eye(d) == 1) = 1;
r = 0.03;
del = 0.7578;
T = 1;
w = 0.25*ones(1,4);

K = [0.5:0.025:1.5]'*(w*S'); 

filename1 = 'basket_K_change.xlsx';
data = readtable(filename1);
LB_data = table2array(data(:,6));

bas1_se = zeros(size(K,1),size(del,2)); 
bas1_price = zeros(size(K,1),size(del,2));
bas1_se_LB = zeros(size(K,1),size(del,2)); %每一行是一个篮子在不同K下的se
bas1_price_LB = zeros(size(K,1),size(del,2));
rate = zeros(size(K,1),1);

phibt = @(u) 2*exp(0.5*u^2*T)* normcdf(sqrt(T)*u, 0, 1);

sigi = sig(1);
kappa = sqrt(2*T/pi);
vecb = w.*S*exp(r*T)/( w*S'*exp(r*T));
cov_mat = diag(sig)*corr*diag(sig);
for l = 1:size(del,2)
    delta = [repmat(del(l),1,d)]';
    deltai = delta(1);
    D_delta = sqrt( eye(d) - diag(delta)^2 ); 
    ncorr_mat = [corr, zeros(1,d)';zeros(1,d),1];
    nmu = zeros(1,d+1);
    nrv1 = randn(n, d+1);  % n x d 矩阵，每行是一个 d 维样本

    L = chol(ncorr_mat, 'lower'); 
    nrv1 = nmu + nrv1 * L';  
    rv1 = nrv1(:,1:d);
    rv2 = nrv1(:,d+1);

    nrv12 = randn(n2, d+1);
    nrv12 = nmu + nrv12 * L'; 
    rv12 = nrv12(:,1:d);
    rv22 = nrv12(:,d+1);

    X = ( D_delta*rv1' + delta*abs(rv2') )' ;
    cl = log( 2* normcdf(sig'.*delta*sqrt(T)) ); 
    coef = log(S') + (r-0.5*sig'.*sig')*T - cl; 
    ST =  exp(coef + repmat(sig' ,1,n).*X'*sqrt(T)); 
    X2 = ( D_delta*rv12' + delta*abs(rv22') )' ;
    ST2 =  exp(coef + repmat(sig' ,1,n2).*X2'*sqrt(T)); 
    LA =  prod(ST.^(vecb'), 1)*prod(w.^(vecb), 2); 
    LA2 =  prod(ST2.^(vecb'), 1)*prod(w.^(vecb), 2); 

    li = log(2* normcdf(sigi*deltai*sqrt(T),  0,  1) ) ;
    term5 = w.*S*exp( (r - 0.5*sigi^2)*T - li );

    for i = 1:length(K)
        LB = LB_data(i);
        L1 = K(i)/( w*S'*exp(r*T))*prod(term5.^vecb,2)*exp(0.5*T*(1-deltai^2)*( sum(sum(vecb'*vecb.*cov_mat)))) * phibt( vecb*sig'*deltai );
            
        contro_LB = (w*ST2- K(i)).*( LA2 >= (L1) ) *exp(-r*T);
        payoff2 = exp(-r*T).*max(0, w*ST2- K(i));
        covLB = cov(contro_LB, payoff2);
        rate(i) = -covLB(1,2)/var(contro_LB);

        payoff_LB = exp(-r*T).*( max(0, w*ST- K(i)) + rate(i)*(w*ST- K(i)).*( LA >= (L1) ) ) - rate(i)*LB; 
        [cp_LB, cv_LB, ci_LB] = normfit(payoff_LB);

        bas1_price_LB(i,l) = cp_LB; 
        bas1_se_LB(i,l) = (ci_LB(2)-ci_LB(1)); 

        payoff = exp(-r*T).*max(0, w*ST- K(i)); 
        [cp, cv, ci] = normfit(payoff);
        bas1_price(i,l) = cp; 
        bas1_se(i,l) = (ci(2)-ci(1)); 
       
    end
end

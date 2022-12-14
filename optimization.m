function ml=optimization(x_0)

load Glasgow
y=[Rentgrowth PE]';
pd= PE;
dg=Rentgrowth;
[K,N]=size(y);
mean_pd=mean(pd);
rho=exp(mean_pd)/(1+exp(mean_pd));
kappa=log(1+exp(mean_pd)) - rho*mean_pd;
 
v=x_0;

%parameters to be estimated
gamma0=v(1);
delta0=v(2);
gamma1=v(3);
delta1=v(4);
sigma_g=abs(v(6));
sigma_mu=abs(v(7));
sigma_d=abs(v(5));
rho_g_mu=v(8);
rho_d_mu=v(9);
rho_g_d=0;
            
%affine coefficients
A=kappa/(1-rho) + (gamma0-delta0)/(1-rho);
B1=1/(1-rho*delta1);
B2=1/(1-rho*gamma1);

Q=sigma_g^2;
R=[sigma_d^2, B2*rho_g_d*sigma_g*sigma_d - B1*rho_d_mu*sigma_mu*sigma_d; ...
    B2*rho_g_d*sigma_g*sigma_d - B1*rho_d_mu*sigma_mu*sigma_d, ...
    B2^2*sigma_g^2 + B1^2*sigma_mu^2 - 2*B1*B2*rho_g_mu*sigma_mu*sigma_g];
S=[rho_g_mu*sigma_g*sigma_mu, B2*sigma_g^2 - B1*rho_g_mu*sigma_g*sigma_mu];


F=gamma1;
G=[1; B2*(gamma1-delta1)];

xpred=0;
omega=sigma_g^2/(1-gamma1^2);
prederrsq=0;
sum_log_delta=0;
ypred=zeros(K,N);
mu=zeros(1,N);

%lap 1
   mu(1)=delta0 + (A+B2*xpred - pd(1))/B1; %%
    delta=G*omega*G'+R;
    
    Theta=F*omega*G' + S;
    
    ypred(:,1)=[gamma0; (1-delta1)*A+delta1*mean_pd] + G*xpred;
    
    inv_delta=inv(delta);

    xpred=F*xpred + Theta*inv_delta*(y(:,1)-ypred(:,1));

    omega=F*omega*F' + Q - Theta*inv_delta*Theta';
    
    prederrsq=prederrsq + (y(:,1)-ypred(:,1))'*inv_delta*(y(:,1)-ypred(:,1));
    sum_log_delta=sum_log_delta + log(det(delta));
    
    
for t=2:N
    
     mu(t)=delta0 + (A+B2*xpred - pd(t))/B1; %%
    delta=G*omega*G'+R;
    
    Theta=F*omega*G' + S;
    
    ypred(:,t)=[gamma0; (1-delta1)*A+delta1*pd(t-1)] + G*xpred;

    inv_delta=inv(delta);

    xpred=F*xpred + Theta*inv_delta*(y(:,t) - ypred(:,t));

    omega=F*omega*F' + Q - Theta*inv_delta*Theta';
    
    prederrsq=prederrsq + (y(:,t)-ypred(:,t))'*inv_delta*(y(:,t)-ypred(:,t));
    sum_log_delta=sum_log_delta + log(det(delta));
    

    
end
    
ml=0.5*(N*K*log(2*pi)+sum_log_delta+prederrsq);
end


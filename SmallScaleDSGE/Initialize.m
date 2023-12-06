
clc
clear all

BETA=0.99;
nu=(1/11);
phi=50;
tau=2;
psi1=1.5;
psi2=0.05;
piss=1
rhog=0.9;
rhoz=0.8;
rhoR=0.85;

save PARAMSTORE  BETA nu phi piss tau psi1 psi2 rhog rhoz rhoR

dynare SmallScale noclearall nolog

omega_ss=oo_.dr.ys;
omega_ss_delta=oo_.dr.ghs2;
F1=oo_.dr.ghx;
F2=oo_.dr.ghu;
F11=oo_.dr.ghxx;
F12=oo_.dr.ghxu;
F22=oo_.dr.ghuu;



M_.endo_names
ordered = M_.endo_names(oo_.dr.order_var,:);

Phi=F1;
R=F2;
Phi11=F11;
Phi12=F12;
R00=F22;
%e(t) ~ iid N(0,Se) 
g=1;
z=2;
Se=zeros(M_.exo_nbr,M_.exo_nbr);
Se(a,a)=sigma_a;
Se(e,e)=sigma_e;



BETA*(exp(c(+1)/exp(c)))^(-tau)*exp(R)/exp(pi(+1))=1;
1=phi*(exp(pi)-1)*((1-0.5*nu)*exp(pi)+0.5*nu)-phi*BETA*(exp(c(+1)/exp(c)))^(-tau)*(exp(y(+1)/exp(y)))*(exp(pi(+1))-1)*exp(pi(+1))+(1-exp(c)^tau)/nu;

R=rhoR*R(-1)+(1-rhoR)*psi1*pi+(1-rhoR)*psi2*(y-g)+u_r;
//g=rhog*g(-1)+u_g;
//z=rhoz*z(-1)+u_z;*/

/*BETA*(exp(c(+1)/exp(c)))^(-tau)*exp(R)/exp(pi(+1))=1;
1=phi*(exp(pi)-1)*((1-0.5*nu)*exp(pi)+0.5*nu)-phi*BETA*(exp(c(+1)/exp(c)))^(-tau)*(exp(y(+1)/exp(y)))*(exp(pi(+1))-1)*exp(pi(+1))+(1-exp(c)^tau)/nu;

R=rhoR*R(-1)+(1-rhoR)*psi1*pi+(1-rhoR)*psi2*(y-g)+u_r;
//g=rhog*g(-1)+u_g;
//z=rhoz*z(-1)+u_z;*/
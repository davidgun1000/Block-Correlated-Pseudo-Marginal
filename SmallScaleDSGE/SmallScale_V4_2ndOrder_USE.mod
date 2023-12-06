
var 
y (long_name='GDP')
x 
c (long_name='Consumption')
R (long_name='Nominal')
real (long_name='Real')
infl (long_name='Inflation')
g (long_name='Fiscal')
z  (long_name='Deamnd')
//x1
;
varexo u_r u_g u_z ;
//  ;

parameters   BETA nu phi piss tau psi1 psi2 rhog  rhoR gss yss gammaQ gamma zss zeta piA  rhoz psi11 psi12;
// r 


load PARAMSTORE;
set_param_value('BETA',BETA);
set_param_value('nu',nu);
set_param_value('tau',tau);
set_param_value('rhoR',rhoR);
set_param_value('psi1',psi1);
set_param_value('psi2',psi2);
set_param_value('phi',phi);
set_param_value('rhog',rhog);
set_param_value('rhoz',rhoz);
set_param_value('zeta',zeta);
set_param_value('gammaQ',gammaQ);
set_param_value('piA',piA);
gss=log(1/(1-zeta));
piss=log(1+piA/400);
yss=log(((1-nu)^(1/tau))/(1-zeta));
gamma=1+gammaQ/100;

//zss=-log(gamma);
zss=log(gamma);

r=log(gamma)-log(BETA);
psi11=(1-rhoR)*psi1;
psi12=(1-rhoR)*psi2;

//gamma=1.02;

model;
//Equation 22 in the write up
//BETA*(exp(c(+1))/exp(c))^(-tau)*(exp(z(+1)))*(exp(R))/exp(infl(+1))=1;
BETA*(exp(c(+1))/exp(c))^(-tau)*(exp(R))/(exp(infl(+1))*exp(z(+1)))=1;

//Equation 23 in the write up
phi*(exp(infl)-exp(piss))*((1-0.5/nu)*exp(infl)+0.5*exp(piss)/nu)-
phi*BETA*(exp(infl(+1))-exp(piss))*exp(infl(+1))*(exp(c(+1))/exp(c))^(-tau)*(exp(y(+1))/exp(y))
+1/nu*(1-(exp(c))^(tau))=1;

//Equation 24 in the write up
exp(R)=exp(real)^(1-rhoR)*exp(R(-1))^rhoR*exp(u_r);
//Equation 25 in the write up
//exp(real)=((gamma)/BETA)*exp(piss)*(exp(infl)/exp(piss))^psi1*(exp(y)/exp(yss))^psi2;
exp(real)=((gamma)/BETA)*exp(piss)*(exp(infl)/exp(piss))^psi1*(exp(y)/(((1-nu)^(1/tau))*exp(g)))^psi2;

/*
//Taylor Rule form --equation 30 in write up
R=(1-rhoR)*(r+piss)+rhoR*R(-1)+psi11*(infl-piss)+psi12*(y-yss)+u_r;*/
//log(exp(R)/exp(r+piss))  = rhoR*log(exp(R(-1))/exp(r+piss))+(1-rhoR)*(BETTAPAI*log(exp(pai)/PAI)+BETTAY*log(exp(y)/Y)
//log(exp(r)/R)  = RHOR*log(exp(r(-1))/R)+(1-RHOR)*(BETTAPAI*log(exp(pai)/PAI)+BETTAY*log(exp(y)/Y)
              

//Equation 26 in the write up
exp(c)/exp(y)=1/exp(g)-(phi*0.5*exp(piss)^2)*(exp(infl)/exp(piss)-1)^2;
//Equation 28 in the write up
g=(1-rhog)*gss+rhog*g(-1)+u_g;
//Equation 29 in the write up
z=(1-rhoz)*zss+rhoz*z(-1)+u_z;
//z=zss;
//exp(x)=exp(R);
exp(x)=exp(y(-1));
end;

initval;
//g=log(1/zeta);
//z=zss;
c=log((1-nu)^(1/tau));
y=yss;
/*R=log(gamma)-log(BETA)+piss;*/

real=log((gamma)/BETA);
R=log((gamma)*exp(piss)/BETA);
infl=log(1+piA/400);
end;


shocks;
var u_r; stderr sigma_r;
var u_g; stderr sigma_g;
var u_z; stderr sigma_z;
end;


options_.noprint=1;
steady;
check;





stoch_simul(order = 2,pruning, irf=0,  nomoments); //Second order

//stoch_simul(order = 1 irf=0,  nomoments); //First order




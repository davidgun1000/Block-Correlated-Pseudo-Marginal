
var 
lambdatilde
ct
Vct
mct
pit
pittilde
pitw
pittildew
wt
nt
Vit
it
qt
Ayt
rtk
ut
aut
ktbar
Rt
mut
etat
Gzt
gt
xtg
yt
yt_lag
ct_lag
it_lag
zt
//Flexible price economy variables

Vctf 
ctf 
lambdatildef 
mctf 
wtf 
ntf 
ytf 
qtf 
itf 
Vitf 
autf 
utf 
rtkf 
Rtf
Aytf 
ktbarf
pitf 
pittildef 
pitwf 
pittildewf

;

varexo epsilongt epsilonmut epsilonetat epsilonzt epsilonRt;

parameters BETA_EST BETA alpha DELTA epsilonp epsilonw gamma sigmaL psiI psip psiw rhoR gammapi
gammax gammag a aest aw awest rhog rhoeta rhomu g pibar Nss c_y numerator denominator psiL Rbar
sigma_g sigma_mu sigma_eta sigma_z sigma_r sigma_gest sigma_muest sigma_etaest sigma_zest sigma_rest
eta mu rk sigmaa Gzest Gz pibar_est 
psipf psiwf pibarf


;

load PARAMSTORE;
set_param_value('BETA_EST',BETA_EST);
set_param_value('Gzest',Gzest);
set_param_value('pibar_est',pibar_est);

set_param_value('alpha',alpha);
set_param_value('DELTA',DELTA);
set_param_value('epsilonp',epsilonp);
set_param_value('epsilonw',epsilonw);
set_param_value('gamma',gamma);
set_param_value('sigmaL',sigmaL);
set_param_value('psiI',psiI);
set_param_value('psip',psip);
set_param_value('psiw',psiw);
set_param_value('rhoR',rhoR);
set_param_value('gammapi',gammapi);
set_param_value('gammax',gammax);
set_param_value('gammag',gammag);
//set_param_value('Nss',Nss);
set_param_value('sigmaa',sigmaa);
set_param_value('psiL',psiL);



//Indexation
set_param_value('aest',aest);
set_param_value('awest',awest);

//AR(1) processes
set_param_value('rhog',rhog);
set_param_value('rhoeta',rhoeta);
set_param_value('rhomu',rhomu);
set_param_value('sigma_gest',sigma_gest);
set_param_value('sigma_muest',sigma_muest);
set_param_value('sigma_etaest',sigma_etaest);
set_param_value('sigma_zest',sigma_zest);
set_param_value('sigma_rest',sigma_rest);

//Shocks and steady state
set_param_value('g',g);
set_param_value('eta',eta);
set_param_value('mu',mu);

pibar=(pibar_est/100+1);
Gz=exp(Gzest/100);
BETA=(1+BETA_EST/100)^-1;

a=1-aest;
aw=1-awest;

sigma_g=sigma_gest/100;
sigma_mu=sigma_muest/100;
sigma_eta=sigma_etaest/100;
sigma_z=sigma_zest/100;
sigma_r=sigma_rest/100;

//psiL=1;//Following Gustt et al
mc=(epsilonp-1)/epsilonp;
c_y=1/g-(1-Gz^(-1)*(1-DELTA))*mc*alpha*Gz/(BETA^(-1)*Gz-1+DELTA);
numerator=((epsilonw-1)/epsilonw)*mc*(1-alpha)*(1-BETA*gamma/Gz);
denominator=(1-gamma/Gz)*c_y;
Nss=(numerator/denominator)^(1/(1+sigmaL));
//denominator=Nss^(1+sigmaL)*(1-gamma/Gz)*c_y;
//psiL=numerator/denominator;
Rbar=BETA^(-1)*Gz*pibar;
rk=BETA^(-1)*Gz-1+DELTA;



psipf=0;
psiwf=0;
pibarf=pibar;

model;
//Recursive expression for marginal utility
Vct=(1/Gzt(+1))*(ct(+1)-gamma*ct/Gzt(+1))^(-1);

//Marginal utility with internal habits
lambdatilde=(ct-gamma*ct(-1)/Gzt)^(-1)-BETA*gamma*Vct;

//RMC
mct=wt*nt/((1-alpha)*yt);
//mct=wt*nt^(alpha)*(Gzt*ktbar)^(-alpha)/(1-alpha);

//Price Phillips curve implied by rotemberg
//(pit/pittilde(-1)-1)*pit/pittilde(-1)=BETA*(lambdatilde(+1)/lambdatilde)*(yt(+1)/yt)*(pit(+1)/pittilde-1)*pit(+1)/pittilde+(epsilonp/psip)*(mct-(epsilonp-1)/epsilonp);
psip*(pit/pittilde(-1)-1)*pit/pittilde(-1)=psip*BETA*(lambdatilde(+1)/lambdatilde)*(yt(+1)/yt)*(pit(+1)/pittilde-1)*pit(+1)/pittilde+(epsilonp/1)*(mct-(epsilonp-1)/epsilonp);

//Wage Phillips curve implied by rotemberg 5
//(pitw/pittildew-1)*pitw/pittildew=BETA*(pitw(+1)/pittildew(+1)-1)*pitw(+1)/pittildew(+1)+lambdatilde*nt*(epsilonw/psiw)*(psiL*nt^sigmaL/lambdatilde-wt*(epsilonw-1)/epsilonw);
psiw*(pitw/pittildew-1)*pitw/pittildew=psiw*BETA*(pitw(+1)/pittildew(+1)-1)*pitw(+1)/pittildew(+1)+lambdatilde*nt*(epsilonw/1)*(psiL*nt^sigmaL/lambdatilde-wt*(epsilonw-1)/epsilonw);

//Optimal condition for investment
1=qt*exp(mut)*(1-psiI*(exp(epsilonzt)*it/it(-1)-1)*exp(epsilonzt)*it/it(-1))+Vit;
Vit=BETA*psiI*qt(+1)*exp(mut(+1))*(lambdatilde(+1)/lambdatilde)*(exp(epsilonzt(+1))*it(+1)/it-1)*(it(+1)^2/it^2)*exp(epsilonzt(+1))-qt*exp(mut)*psiI*0.5*(exp(epsilonzt)*it/it(-1)-1)^2;

//Optimal condition for capital

aut=rk/sigmaa*(exp(sigmaa*(ut-1))-1);

qt=BETA*(lambdatilde(+1)/lambdatilde/Gzt(+1))*(qt(+1)*(1-DELTA)+(rtk(+1))*ut(+1)-aut(+1));
//1=BETA*(lambdatilde(+1)/lambdatilde/Gzt(+1))*(1+rt(+1));
//rt(+1)=(qt(+1)*(1-DELTA)+(rtk(+1))*ut(+1)-aut(+1))/qt-1;

rtk=rk*(exp(sigmaa*(ut-1)));

//rtk=alpha/(1-alpha)*Gzt*wt*nt/(ut*ktbar); 
rtk=alpha/(1-alpha)*Gzt*wt*nt/(ut*ktbar(-1)); 

wt=pitw*wt(-1)/(pit*Gzt);

//pittilde(-1)=pibar^a*(pit(-1))^(1-a);
pittilde=pibar^a*(pit)^(1-a);

pittildew=Gz*(pibar^aw)*(exp(epsilonzt)*pit(-1))^(1-aw);

//Euler equation
lambdatilde=BETA*exp(etat)*Rt*(lambdatilde(+1)/(Gzt(+1)*pit(+1)));

//Market clearing
//yt=Ayt^(-1)*(ct+it+aut*ktbar/Gzt);
yt=Ayt^(-1)*(ct+it+aut*ktbar(-1)/Gzt);

Ayt=1/(gt)-psip*0.5*(pit/pittilde(-1)-1)^2;

//Labor choice
//nt=yt^(1/(1-alpha))*(ut*ktbar/Gzt)^(alpha/(alpha-1));
//yt=(ut*ktbar/Gzt)^alpha*nt^(1-alpha);
yt=(ut*ktbar(-1)/Gzt)^alpha*nt^(1-alpha);

//Law of motion of capital
//ktbar(+1)=(1-DELTA)*ktbar*Gzt^(-1)+exp(mut)*(1-psiI*0.5*(it*exp(epsilonzt)/it(-1)-1)^2)*it;
ktbar=(1-DELTA)*ktbar(-1)*Gzt^(-1)+exp(mut)*(1-psiI*0.5*(it*exp(epsilonzt)/it(-1)-1)^2)*it;

//Monetary policy Taylor rule
log(Rt/Rbar)=rhoR*log(Rt(-1)/Rbar)+(1-rhoR)*(gammapi*log(pit/pibar)+gammax*log(yt/ytf)+gammag*log(yt*exp(epsilonzt)/yt(-1))) +epsilonRt;


//Natural output
(xtg)=alpha*log(ut)+(1-alpha)*log(nt/Nss);
//exogenous processes
exp(gt)=exp(g)^(1-rhog)*exp(gt(-1))^(rhog)*exp(epsilongt);
//gt=(1-rhog)*g+rhog*gt(-1)+epsilongt;
exp(mut)=exp(mu)^(1-rhomu)*exp(mut(-1))^(rhomu)*exp(epsilonmut);
//mut=rhomu*mut(-1)+epsilonmut;

//exp(etat)=exp(eta)^(1-rhoeta)*exp(etat(-1))^(rhoeta)*exp(epsilonetat);
etat=rhoeta*etat(-1)+epsilonetat;

Gzt=Gz*exp(epsilonzt);
yt_lag=yt(-1);
ct_lag=ct(-1);
it_lag=it(-1);
zt=epsilonzt;


//Flexible price economy

Vctf=(1/Gzt(+1))*(ctf(+1)-gamma*ctf/Gzt(+1))^(-1);

//Marginal utility with internal habits
lambdatildef=(ctf-gamma*ctf(-1)/Gzt)^(-1)-BETA*gamma*Vctf;

//RMC
mctf=wtf*ntf/((1-alpha)*ytf);

//Price Phillips curve implied by rotemberg
//0=(mctf-(epsilonp-1)/epsilonp);
psipf*(pitf/pittildef(-1)-1)*pitf/pittildef(-1)=psipf*BETA*(lambdatildef(+1)/lambdatildef)*(ytf(+1)/ytf)*(pitf(+1)/pittildef-1)*pitf(+1)/pittildef+(epsilonp/1)*(mctf-(epsilonp-1)/epsilonp);


//Wage Phillips curve implied by rotemberg 5
//0=(psiLf*ntf^sigmaL/lambdatildef-wtf*(epsilonw-1)/epsilonw);
psiwf*(pitwf/pittildewf-1)*pitwf/pittildewf=psiwf*BETA*(pitwf(+1)/pittildewf(+1)-1)*pitwf(+1)/pittildewf(+1)+lambdatildef*ntf*(epsilonw/1)*(psiL*ntf^sigmaL/lambdatildef-wtf*(epsilonw-1)/epsilonw);

wtf=pitwf*wtf(-1)/(pitf*Gzt);

pittildef=pibarf^a*(pitf)^(1-a);

pittildewf=Gz*(pibarf^aw)*(exp(epsilonzt)*pitf(-1))^(1-aw);

rtkf=alpha/(1-alpha)*Gzt*wtf*ntf/(utf*ktbarf(-1)); 


//Optimal condition for investment
1=qtf*exp(mut)*(1-psiI*(exp(epsilonzt)*itf/itf(-1)-1)*exp(epsilonzt)*itf/itf(-1))+Vitf;

Vitf=BETA*psiI*qtf(+1)*exp(mut(+1))*(lambdatildef(+1)/lambdatildef)*(exp(epsilonzt(+1))*itf(+1)/itf-1)*(itf(+1)^2/itf^2)*exp(epsilonzt(+1))-qtf*exp(mut)*psiI*0.5*(exp(epsilonzt)*itf/itf(-1)-1)^2;

//Optimal condition for capital

autf=rk/sigmaa*(exp(sigmaa*(utf-1))-1);

qtf=BETA*(lambdatildef(+1)/lambdatildef/Gzt(+1))*(qtf(+1)*(1-DELTA)+(rtkf(+1))*utf(+1)-autf(+1));

rtkf=rk*(exp(sigmaa*(utf-1)));

//Euler equation
lambdatildef=BETA*exp(etat)*Rtf*(lambdatildef(+1)/(Gzt(+1)*pitf(+1)));

//Market clearing
ytf=Aytf^(-1)*(ctf+itf+autf*ktbarf(-1)/Gzt);

Aytf=1/(gt);

ytf=(utf*ktbarf(-1)/Gzt)^alpha*ntf^(1-alpha);

ktbarf=(1-DELTA)*ktbarf(-1)*Gzt^(-1)+exp(mut)*(1-psiI*0.5*(itf*exp(epsilonzt)/itf(-1)-1)^2)*itf;

//START HERE WITH THE TALOR RULE
log(Rtf/Rbar)=rhoR*log(Rtf(-1)/Rbar)+(1-rhoR)*(gammapi*log(pitf/pibar)+gammax*log(ytf/ytf)+gammag*log(ytf*exp(epsilonzt)/ytf(-1))) +epsilonRt;

end;

 


initval;
Gzt=Gz;
pit=pibar;
pittilde=pibar;
mct=(epsilonp-1)/epsilonp;
pittildew=Gzt*pibar;
pitw=Gzt*pibar;
ut=1;
aut=0;
//Rbar=BETA^(-1)*(Gzt*pibar);
qt=1;
rtk=BETA^(-1)*Gzt-1+DELTA;
nt=Nss;
ktbar=Nss*Gzt*((alpha/(1-alpha))*mc/(BETA^(-1)*Gzt-1+DELTA))^(1/(1-alpha));
//ktbar=alpha/(1-alpha)*Gzt*wt*nt/(BETA^(-1)*Gzt-1+DELTA);
yt=((Gzt^(-1)*ktbar)^alpha)*nt^(1-alpha);
wt=(mc/nt)*(1-alpha)*yt;
it=ktbar*(1-(1-DELTA)/Gzt);
ct=yt/g-it/yt;
//ct=c_y*yt;
Ayt=1/g;
Rt=BETA^(-1)*(Gzt*pibar);
Vct=Gzt^(-1)*(ct-gamma*ct/Gzt)^(-1);
Vit=0;
lambdatilde=(1-Gzt^(-1)*gamma)^(-1)*(1-BETA*Gzt^(-1)*gamma)/ct;
gt=g;
xtg=0;
mut=0;
yt_lag=yt;
ct_lag=ct;
it_lag=it;


ctf=ct;
lambdatildef=lambdatilde;
mctf=mct;
wtf=wt;
ntf=Nss;
ytf=yt;
qtf=qt;
itf=it;
Vitf=Vit;
autf=aut;
utf=ut;
rtkf=rtk;
Rtf=Rt;
Aytf=Ayt;
ktbarf=ktbar;
pitf=pit;
pittildef=pittilde;
pitwf=pitwf;
pittildewf=pittildew;
end;

shocks;
var epsilongt; stderr sigma_g;
var epsilonmut; stderr sigma_mu;
var epsilonetat; stderr sigma_eta;
var epsilonzt; stderr sigma_z;
var epsilonRt;stderr sigma_r;

end;

options_.noprint=1;
steady(solve_algo=1,maxit=500000);      

[eigenvalues_,result,info] = check(M_, options_, oo_);


if result==1
stoch_simul(order = 2,pruning,  irf=0, nomoments); 
estrun=result;
else
estrun=result;
display('no solution')
return
end


//stoch_simul(order = 1 irf=0,  nomoments); //First order



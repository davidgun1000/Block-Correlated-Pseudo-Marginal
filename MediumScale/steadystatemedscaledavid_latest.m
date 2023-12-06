
clc
clear all
addpath 'C:/dynare/4.4.3/matlab'
% addpath 'C:/dynare/4.5.7/matlab'

%addpath 'C:\Users\prati\My Drive\dynare\4.4.3\matlab'
%addpath 'C:\Users\z3526523\Dropbox\NonlinearEstimation\Revision\MediumScale\Debug\SentPratitiError'

%Calibrated parameters
DELTA=0.025;%Depreciation rate of capital
g=1.25;%share of government spending 20% i.e. 1-1/g=0.2 in steady state
epsilonp=6;%steady state net markup: epsilonp/(epsilonp-1)-1=0.2
epsilonw=6;%steady state net markup: epsilonw/(epsilonw-1)-1=0.2
rhoeta=0.85;%Persistence of liquidity shock 
psiL=1;

%Estimated Parameters
% BETA=0.99; %Need to transform the parameter expression
% BETA=(1+BETA_EST/100)^-1;

%%USE THESE INITIAL VALUES
BETA_EST=100*(0.995^-1-1);
% Gz=exp(Gzest/100);
Gzest=100*log(1.0050);
% pibar=1.02;
%pibar=(pibar_est/100+1)
pibar_est=100*(1.0061-1);
% pibar_est=100*(1.0126-1);

alpha=1/3;
%Endogenous propagation
gamma=0.6;
sigmaa=5;
psiI=5;
psip=80;
psiw=3000;
sigmaL=4; 
%Price and wage indexation
% a=1-aest;
aest=0.5;
%aw=1-aest;
awest=0.5;


%Monetary Policy
rhoR=0.9;
gammapi=2;
gammax=0.14;
gammag=0.7;

%AR(1) coefficients for shocks
rhog=0.5;
rhomu=0.6;
%standard deviation for shocks
sigma_gest=100*0.002;
sigma_muest=100*0.002;
sigma_etaest=100*0.002;
sigma_zest=100*0.002;
sigma_rest=100*0.002;
eta=0;
mu=0;

% Nss=1;





save PARAMSTORE BETA_EST alpha DELTA epsilonp epsilonw gamma sigmaL psiI psip psiw rhoR gammapi... 
gammax gammag aest awest rhog rhoeta rhomu g Gzest pibar_est   ...
sigma_gest sigma_muest sigma_etaest sigma_zest sigma_rest eta mu sigmaa psiL;


%states in the model
%ktbar yt ct it pit Rt wt eta gt yt(-1) ct(-1) it(-1)

% dynare mediumscalenewSW noclearall
dynare medscalenew noclearall

if estrun==1

f_11=[oo_.dr.ghx oo_.dr.ghu];
order=2;

optPruning=[];

[omega_ss,F1,F2,omega_ss_delta,F11,F22,F12,~,~,~,~,~,state_var_positions,control_var_postion] = statespaceform(optPruning,oo_,M_,f_11,order);
Phi=F1;
R=F2;
%Extra Coefficients for seecond order solution
Phi11=F11;
Phi12=F12;
R00=F22;

else 
    return
end

M_.endo_names

ordered = M_.endo_names(oo_.dr.order_var,:);
ordered=cellstr(ordered);

M_.endo_names1=cellstr(M_.endo_names);
%Number of state vriables
M_.endo_names1(oo_.dr.state_var)


obs_ind=[find(strcmp(ordered,'yt')),find(strcmp(ordered,'yt_lag')),find(strcmp(ordered,'ct')),find(strcmp(ordered,'ct_lag')),find(strcmp(ordered,'it')),find(strcmp(ordered,'it_lag')),find(strcmp(ordered,'pit')),find(strcmp(ordered,'Rt')),...
    find(strcmp(ordered,'zt'))];

%Observables -->
%1-->gdp
eq_y=1;
%2--> conusmption
eq_c=2;
%3--> investment
eq_i=3;
%4--> inflation
eq_pi=4;
%5--> ffr
eq_ffr=5;


ny=eq_ffr;

A         = zeros(ny,1);
A(eq_y,1) = Gzest/100;
A(eq_c,1) = Gzest/100;
A(eq_i,1) = Gzest/100;
A(eq_pi,1) = 0;
A(eq_ffr,1) = 0;

nstate = size(Phi,1); 


B = zeros(ny,nstate);

B(eq_y,obs_ind(1)) =  1;
B(eq_y,obs_ind(2)) =  -1; 
B(eq_y,obs_ind(9)) = 1;

B(eq_c,obs_ind(3)) =  1;
B(eq_c,obs_ind(4)) =  -1; 
B(eq_c,obs_ind(9)) = 1;

B(eq_i,obs_ind(5)) =  1;
B(eq_i,obs_ind(6)) =  -1; 
B(eq_i,obs_ind(9)) = 1;

B(eq_pi,obs_ind(7)) =  1;
B(eq_ffr,obs_ind(8)) = 1; 

%Measurement error (ranges from 10% to 25% of total variance observed in
%the data
% data=xlsread('C:\Users\z3526523\Dropbox\NonlinearEstimation\Revision\MediumScale\Data\datause','Sheet1');
%data=xlsread('C:\Users\z3526523\Dropbox\NonlinearEstimation\Revision\MediumScale\Data\datause','Sheet2');
data=xlsread('datause','Sheet2');

%USE THE WU XIA shadow rate in sheet 2.

cc=x2mdate(data(:,1));
cc=datevec(cc);
if cc(1,2)==1
    date_init=cc(1,1)+(1-1)*0.25;
elseif cc(1,2)==4
    date_init=cc(1,1)+(2-1)*0.25;
elseif cc(1,2)==7
    date_init=cc(1,1)+(3-1)*0.25;
elseif cc(1,2)==10
    date_init=cc(1,1)+(4-1)*0.25;
end
date(1)=date_init;
for i=1:length(data)-1
    date(i+1,1)=date(i)+0.25;
end

%cutting of sample before zero lower bound on nominal interest rates
%Sample 1
% endsample=find(date(:,1)==[2007.75000000000]);
%Sample 2
endsample=find(date(:,1)==[2019.75000000000]);

datasub=data(1:endsample,2:end);

m_E=0.1; %MEASUREMENT ERROR OF 10% (Gust et al also try with 25% of total variance observed in
%the data; m_E=0.25)

Sigma_u=m_E*diag(var(datasub));
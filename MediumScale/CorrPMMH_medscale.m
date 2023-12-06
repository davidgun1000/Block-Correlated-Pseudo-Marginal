%simulate med scale

clear
clc
close all
delete *.asv
tic
%addpath 'C:\dynare\4.6.0\matlab'
%addpath 'C:\Users\dgunawan\Dropbox\Gunawan_Kohn\DavidCodev1\RBC_sim'
%addpath 'C:\dynare\4.4.3\matlab'
addpath 'C:\dynare\4.5.7\matlab'
%addpath 'C:\dynare\5.1\matlab'

%addpath 'D:\Dropbox\Gunawan_Kohn\DavidCodev1\RBC_sim'
%addpath '/Applications/Dynare/4.5.7/matlab'

%For NIASRA
addpath('/home/davidg/dynare/lib/dynare/matlab')

%For KATANA
%addpath '/apps/dynare/4.4.3/lib/dynare/matlab'


%for GADI
addpath('/home/561/dg2271/dynare457/lib/dynare/matlab')
addpath('/home/561/dg2271/matio1.5.17') 
%parpool(48)
load('cov_medscale.mat');
load('medscale_data_used.mat');
yt=dataused;

DELTA=0.025;%Depreciation rate of capital
g=1.25;%share of government spending 20% i.e. 1-1/g=0.2 in steady state
epsilonp=6;%steady state net markup: epsilonp/(epsilonp-1)-1=0.2
epsilonw=6;%steady state net markup: epsilonw/(epsilonw-1)-1=0.2
rhoeta=0.85;%Persistence of liquidity shock 
psiL=1;

scale_psip = 100;
scale_psiw = 1000;
scale_param_shock=100;
param_BETA_EST=0.0610;
param_Gzest=0.5282;
param_pibar_est = 0.5756;
param_alpha=0.0933;
param_gamma=0.7157;
param_sigmaa=5.9849;
param_psiI=3.6219;
param_psip_bar=43.6217/scale_psip;
param_psiw_bar=6.2646e+03/scale_psiw;
param_sigmaL= 4.9567;
param_aest=0.5763;
param_awest=0.6199;
param_rhoR=0.6301;
param_gammapi=1.0321;
param_gammax=0.2664;
param_gammag=0.7433;
param_rhog=0.7738;
param_rhomu=0.8861;
param_sigma_gest=0.0955;
param_sigma_muest=3.3335;
param_sigma_etaest=0.5138;
param_sigma_zest=0.4996;
param_sigma_rest=0.1918;
param_eta=0;
param_mu=0;
%param_Nss=1;
param_sigma_measure1_bar = scale_param_shock*(0.0036);
param_sigma_measure2_bar = scale_param_shock*(0.0038);
param_sigma_measure3_bar = scale_param_shock*(0.0275);
param_sigma_measure4_bar = scale_param_shock*(0.0012);
param_sigma_measure5_bar = scale_param_shock*(0.0035);

BETA_EST=param_BETA_EST;
Gzest=param_Gzest;
pibar_est = param_pibar_est;
alpha=param_alpha;
gamma=param_gamma;
sigmaa=param_sigmaa;
psiI=param_psiI;
psip=scale_psip*param_psip_bar;
psiw=scale_psiw*param_psiw_bar;
sigmaL=param_sigmaL;
aest=param_aest;
awest=param_awest;
rhoR=param_rhoR;
gammapi=param_gammapi;
gammax=param_gammax;
gammag=param_gammag;
rhog=param_rhog;
rhomu=param_rhomu;
sigma_gest=param_sigma_gest;
sigma_muest=param_sigma_muest;
sigma_etaest=param_sigma_etaest;
sigma_zest=param_sigma_zest;
sigma_rest=param_sigma_rest;
eta=param_eta;
mu=param_mu;
%Nss=param_Nss;
sigma_measure1 = param_sigma_measure1_bar/scale_param_shock;
sigma_measure2 = param_sigma_measure2_bar/scale_param_shock;
sigma_measure3 = param_sigma_measure3_bar/scale_param_shock;
sigma_measure4 = param_sigma_measure4_bar/scale_param_shock;
sigma_measure5 = param_sigma_measure5_bar/scale_param_shock;

save PARAMSTORE BETA_EST alpha DELTA epsilonp epsilonw gamma sigmaL psiI psip psiw rhoR gammapi... 
gammax gammag aest awest rhog rhoeta rhomu g Gzest pibar_est   ...
sigma_gest sigma_muest sigma_etaest sigma_zest sigma_rest eta mu sigmaa psiL;


%dynare eig_check noclearall nolog
%states in the model
%ktbar yt ct it pit Rt wt eta gt yt(-1) ct(-1) it(-1)

dynare medscalenew noclearall nolog

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
delta_ss = omega_ss_delta;
S2 = M_.Sigma_e;

nstate = size(Phi,1); 
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

B = B_true;
H=diag([sigma_measure1^2;sigma_measure2^2;sigma_measure3^2;sigma_measure4^2;sigma_measure5^2]);


%prior distribution hyperparameters
prior_rho_mean=0.6;
prior_rho_std=0.2;
prior_gamma_a=0.5;
prior_gamma_b=2;
prior_gammag_1 = 0.4;
prior_gammag_2 = 0.3;
prior_gammapi_1 = 1.7;
prior_gammapi_2 = 0.3;
prior_gammax_1 = 0.4;
prior_gammax_2 = 0.3;

prior_gamma_1 = 0.6;
prior_gamma_2 = 0.1;
prior_sigmaa_1 = 25;
prior_sigmaa_2 = 5;

prior_psip_1 = 1;
prior_psip_2 = 0.25;
prior_psiw_1 = 3;
prior_psiw_2 = 5;
prior_psiI_1 = 16;
prior_psiI_2 = 4;

prior_sigmaL_1 = 5.34;
prior_sigmaL_2 = 2.66;

prior_a_1 = 0.5;
prior_a_2 = 0.15;

prior_aw_1 = 0.5;
prior_aw_2 = 0.15;

prior_BETAEST_1 = 0.625;
prior_BETAEST_2 = 2.5;

prior_Gzest_1 = 0.5;
prior_Gzest_2 = 0.03;

prior_pibarest_1 = 0.62;
prior_pibarest_2 = 0.1;

prior_alpha_1 = 0.3;
prior_alpha_2 = 0.05;
 
S=1; %number of nodes
N=50000;
burn=2000;
adapt=3000;
nit=200000;
nloop=burn+nit+adapt;
ne=size(S2,1);
[~, ns] = size(B);
T=size(yt,2);
num_prop=100;
N_prop=1000;
corr_param=0.9;
num_param=28;
num_mix=1;
corr_index=1;
disturbance_index=1;
w_prior=0.05;
w_rest=1-w_prior;

for s=1:S    
    u_particles(:,1:N,:,s)=randn(ne,N,T);
    u_res(:,1:N,s)=rand(T,N);
end
for s=1:S
    [particles_eps(:,:,:,s),particles_first(:,:,:,s),particles_second(:,:,:,s),w(:,:,s),indx(:,:,s),lik_indv(:,s)]=Mix_smc_SmallScale_2ndorder_bootstrap_corr(A,B,H,Phi,R,S2,Phi11,Phi12,R00,delta_ss,N,yt,u_particles(:,:,:,s),u_res(:,:,s),...
        state_var_positions,corr_index,disturbance_index);
end
lik_PF=sum(lik_indv); % average of the log likelihood

prior_PF=log_truncatednormal(param_rhoR,prior_rho_mean,prior_rho_std,0,1)+...
      log_truncatednormal(param_rhog,prior_rho_mean,prior_rho_std,0,1)+log_truncatednormal(param_rhomu,prior_rho_mean,prior_rho_std,0,1)+...
      log(gampdf(param_sigma_gest,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_muest,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_etaest,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_zest,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_rest,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_measure1_bar,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_measure2_bar,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_measure3_bar,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_measure4_bar,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_measure5_bar,prior_gamma_a,1/prior_gamma_b))+...
      log(normpdf(param_gammag,prior_gammag_1,prior_gammag_2))+log(normpdf(param_gammapi,prior_gammapi_1,prior_gammapi_2))+...
      log(normpdf(param_gammax,prior_gammax_1,prior_gammax_2))+...
      log_truncatednormal(param_gamma,prior_gamma_1,prior_gamma_2,0,1)+log(gampdf(param_sigmaa,prior_sigmaa_1,1/prior_sigmaa_2))+...
      log(normpdf(param_psip_bar,prior_psip_1,prior_psip_2))+log(normpdf(param_psiw_bar,prior_psiw_1,prior_psiw_2))+...
      log(gampdf(param_psiI,prior_psiI_1,1/prior_psiI_2))+log(gampdf(param_sigmaL,prior_sigmaL_1,1/prior_sigmaL_2))+...
      log_truncatednormal(param_aest,prior_a_1,prior_a_2,0,1)+log_truncatednormal(param_awest,prior_aw_1,prior_aw_2,0,1)+...
      log(gampdf(param_BETA_EST,prior_BETAEST_1,1/prior_BETAEST_2))+log(normpdf(param_Gzest,prior_Gzest_1,prior_Gzest_2))+...
      log(normpdf(param_pibar_est,prior_pibarest_1,prior_pibarest_2))+log(normpdf(param_alpha,prior_alpha_1,prior_alpha_2));

      
post_PF=prior_PF+lik_PF;

jac_PF=log((1/param_rhoR)+(1/(1-param_rhoR)))+log((1/param_rhog)+(1/(1-param_rhog)))+log((1/param_rhomu)+(1/(1-param_rhomu)))+...
    log(1/param_sigma_gest)+log(1/param_sigma_muest)+log(1/param_sigma_etaest)+log(1/param_sigma_zest)+log(1/param_sigma_rest)+...
    log(1/param_sigma_measure1_bar)+log(1/param_sigma_measure2_bar)+log(1/param_sigma_measure3_bar)+log(1/param_sigma_measure4_bar)+...
    log(1/param_sigma_measure5_bar)+log((1/param_gamma)+(1/(1-param_gamma)))+log(1/param_sigmaa)+log(1/param_psiI)+log(1/param_sigmaL)+...
    log((1/param_aest)+(1/(1-param_aest)))+log((1/param_awest)+(1/(1-param_awest)))+log(1/param_BETA_EST);

D1=28;
%V1=0.01*eye(D1);
V1=covmat_theta;
scale=0.1;
accept=0;
A_rand=0;
init_time=100;
target_accept=0.20;
r1_PF=0;
r1_CDKF=0;
lik_PF_star=0;
lik_CDKF_star=0;
accept_CDKF=0;

for i=1:nloop
    i
    accept
    %param_sigmaa
    %tic
    rng('shuffle')
    theta=[logit_inverse(param_rhoR),logit_inverse(param_rhog),logit_inverse(param_rhomu),log(param_sigma_gest),...
        log(param_sigma_muest),log(param_sigma_etaest),log(param_sigma_zest),log(param_sigma_rest),...
        log(param_sigma_measure1_bar),log(param_sigma_measure2_bar),log(param_sigma_measure3_bar),log(param_sigma_measure4_bar),...
        log(param_sigma_measure5_bar),param_gammag,param_gammapi,param_gammax,logit_inverse(param_gamma),log(param_sigmaa),...
        param_psip_bar,param_psiw_bar,log(param_psiI),log(param_sigmaL),logit_inverse(param_aest),logit_inverse(param_awest),...
        log(param_BETA_EST),param_Gzest,param_pibar_est,param_alpha]; 
            
    rng('shuffle')
    R1=mvnrnd(theta,scale.*V1);
    
    param_rhoR_star=logitcdf(R1(1,1));
    param_rhog_star=logitcdf(R1(1,2));
    param_rhomu_star=logitcdf(R1(1,3));
    param_sigma_gest_star = exp(R1(1,4));
    param_sigma_muest_star = exp(R1(1,5));
    param_sigma_etaest_star = exp(R1(1,6));
    param_sigma_zest_star = exp(R1(1,7));
    param_sigma_rest_star = exp(R1(1,8));
    param_sigma_measure1_bar_star = exp(R1(1,9));
    param_sigma_measure2_bar_star = exp(R1(1,10));
    param_sigma_measure3_bar_star = exp(R1(1,11));
    param_sigma_measure4_bar_star = exp(R1(1,12));
    param_sigma_measure5_bar_star = exp(R1(1,13));
    param_gammag_star = R1(1,14);
    param_gammapi_star = R1(1,15);
    param_gammax_star = R1(1,16);
    param_gamma_star = logitcdf(R1(1,17));
    param_sigmaa_star = exp(R1(1,18));
    param_psip_bar_star = R1(1,19);
    param_psiw_bar_star = R1(1,20);
    param_psiI_star = exp(R1(1,21));
    param_sigmaL_star = exp(R1(1,22));
    param_aest_star = logitcdf(R1(1,23));
    param_awest_star = logitcdf(R1(1,24));
    param_BETA_EST_star = exp(R1(1,25));
    param_Gzest_star = R1(1,26);
    param_pibar_est_star = R1(1,27);
    param_alpha_star = R1(1,28);
    
    BETA_EST=param_BETA_EST_star;
    Gzest=param_Gzest_star;
    pibar_est = param_pibar_est_star;
    alpha=param_alpha_star;
    gamma=param_gamma_star;
    sigmaa=param_sigmaa_star;
    psiI=param_psiI_star;
    psip=scale_psip*param_psip_bar_star;
    psiw=scale_psiw*param_psiw_bar_star;
    sigmaL=param_sigmaL_star;
    aest=param_aest_star;
    awest=param_awest_star;
    rhoR=param_rhoR_star;
    gammapi=param_gammapi_star;
    gammax=param_gammax_star;
    gammag=param_gammag_star;
    rhog=param_rhog_star;
    rhomu=param_rhomu_star;
    sigma_gest=param_sigma_gest_star;
    sigma_muest=param_sigma_muest_star;
    sigma_etaest=param_sigma_etaest_star;
    sigma_zest=param_sigma_zest_star;
    sigma_rest=param_sigma_rest_star;
    eta=param_eta;
    mu=param_mu;
    %Nss=param_Nss;
    sigma_measure1 = param_sigma_measure1_bar_star/scale_param_shock;
    sigma_measure2 = param_sigma_measure2_bar_star/scale_param_shock;
    sigma_measure3 = param_sigma_measure3_bar_star/scale_param_shock;
    sigma_measure4 = param_sigma_measure4_bar_star/scale_param_shock;
    sigma_measure5 = param_sigma_measure5_bar_star/scale_param_shock;

    save PARAMSTORE BETA_EST alpha DELTA epsilonp epsilonw gamma sigmaL psiI psip psiw rhoR gammapi... 
    gammax gammag aest awest rhog rhoeta rhomu g Gzest pibar_est   ...
    sigma_gest sigma_muest sigma_etaest sigma_zest sigma_rest eta mu sigmaa psiL;

     
    %eig_check
    
    %if size(find(real(abs(oo_.dr.eigval))>1),1)==M_.nsfwrd
    
    medscalenew
        
    if result==1 & param_sigma_gest_star>0.005 & param_sigma_muest_star>0.005 & param_sigma_etaest_star>0.005 & param_sigma_zest_star>0.005 ...
            & param_sigma_rest_star>0.005 & param_rhoR_star>0.001 & param_rhoR_star<0.999 & param_rhog_star>0.001 & param_rhog_star<0.999 ...
            & param_rhomu_star>0.001 & param_rhomu_star<0.999 & param_sigma_measure1_bar_star>0.005 & param_sigma_measure2_bar_star>0.005 ...
            & param_sigma_measure3_bar_star>0.005 & param_sigma_measure4_bar_star>0.005 & param_sigma_measure5_bar_star>0.005
        
    f_11=[oo_.dr.ghx oo_.dr.ghu];
    order=2;

    optPruning=[];

    [omega_ss,F1,F2,omega_ss_delta,F11,F22,F12,~,~,~,~,~,state_var_positions,control_var_postion] = statespaceform(optPruning,oo_,M_,f_11,order);
    
    A_star(eq_y,1) = Gzest/100;
    A_star(eq_c,1) = Gzest/100;
    A_star(eq_i,1) = Gzest/100;
    A_star(eq_pi,1) = 0;
    A_star(eq_ffr,1) = 0;
    B_star=B;
    H_star=diag([sigma_measure1^2;sigma_measure2^2;sigma_measure3^2;sigma_measure4^2;sigma_measure5^2]);
    Phi_star=F1;
    R_star=F2;
    S2_star=M_.Sigma_e;
    Phi11_star=F11;
    Phi12_star=F12;
    R00_star=F22;
    delta_ss_star=omega_ss_delta;
    
    rng('shuffle')
    
    u_particles_star=u_particles;
    u_res_star=u_res;
    cindx=randsample(S,1); % select one block to refresh the random numbers

    u_particles_star(:,:,:,cindx)=corr_param*u_particles(:,:,:,cindx)+sqrt(1-corr_param^2)*randn(ne,N,T);
    u_res_norm=norminv(u_res(:,:,cindx));
    u_res_norm_star=corr_param*u_res_norm+sqrt(1-corr_param^2)*randn(T,N);
    u_res_star(:,:,cindx)=normcdf(u_res_norm_star);


   
               
    for s=1:S
        [particles_eps_star(:,:,:,s),particles_first_star(:,:,:,s),particles_second_star(:,:,:,s),w_star(:,:,s),indx_star(:,:,s),lik_indv_star(:,s)]=Mix_smc_SmallScale_2ndorder_bootstrap_corr(A_star,B_star,H_star,Phi_star,R_star,S2_star,...
            Phi11_star,Phi12_star,R00_star,delta_ss_star,N,yt,u_particles_star(:,:,:,s),u_res_star(:,:,s),...
            state_var_positions,corr_index,disturbance_index);        
    end
    %time(i,1)=toc;
    if sum(isnan(lik_indv_star))>0 | sum(sum(sum(sum(isnan(particles_eps_star),4),3),2))>0 
       
       
       
       
    else
    
    lik_PF_star=sum(lik_indv_star);
    
    
    prior_PF_star=log_truncatednormal(param_rhoR_star,prior_rho_mean,prior_rho_std,0,1)+...
      log_truncatednormal(param_rhog_star,prior_rho_mean,prior_rho_std,0,1)+log_truncatednormal(param_rhomu_star,prior_rho_mean,prior_rho_std,0,1)+...
      log(gampdf(param_sigma_gest_star,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_muest_star,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_etaest_star,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_zest_star,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_rest_star,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_measure1_bar_star,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_measure2_bar_star,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_measure3_bar_star,prior_gamma_a,1/prior_gamma_b))+log(gampdf(param_sigma_measure4_bar_star,prior_gamma_a,1/prior_gamma_b))+...
      log(gampdf(param_sigma_measure5_bar_star,prior_gamma_a,1/prior_gamma_b))+...
      log(normpdf(param_gammag_star,prior_gammag_1,prior_gammag_2))+log(normpdf(param_gammapi_star,prior_gammapi_1,prior_gammapi_2))+...
      log(normpdf(param_gammax_star,prior_gammax_1,prior_gammax_2))+...
      log_truncatednormal(param_gamma_star,prior_gamma_1,prior_gamma_2,0,1)+log(gampdf(param_sigmaa_star,prior_sigmaa_1,1/prior_sigmaa_2))+...
      log(normpdf(param_psip_bar_star,prior_psip_1,prior_psip_2))+log(normpdf(param_psiw_bar_star,prior_psiw_1,prior_psiw_2))+...
      log(gampdf(param_psiI_star,prior_psiI_1,1/prior_psiI_2))+log(gampdf(param_sigmaL_star,prior_sigmaL_1,1/prior_sigmaL_2))+...
      log_truncatednormal(param_aest_star,prior_a_1,prior_a_2,0,1)+log_truncatednormal(param_awest_star,prior_aw_1,prior_aw_2,0,1)+...
      log(gampdf(param_BETA_EST_star,prior_BETAEST_1,1/prior_BETAEST_2))+log(normpdf(param_Gzest_star,prior_Gzest_1,prior_Gzest_2))+...
      log(normpdf(param_pibar_est_star,prior_pibarest_1,prior_pibarest_2))+log(normpdf(param_alpha_star,prior_alpha_1,prior_alpha_2));
    
    jac_PF_star=log((1/param_rhoR_star)+(1/(1-param_rhoR_star)))+log((1/param_rhog_star)+(1/(1-param_rhog_star)))+log((1/param_rhomu_star)+(1/(1-param_rhomu_star)))+...
        log(1/param_sigma_gest_star)+log(1/param_sigma_muest_star)+log(1/param_sigma_etaest_star)+log(1/param_sigma_zest_star)+log(1/param_sigma_rest_star)+...
        log(1/param_sigma_measure1_bar_star)+log(1/param_sigma_measure2_bar_star)+log(1/param_sigma_measure3_bar_star)+log(1/param_sigma_measure4_bar_star)+...
        log(1/param_sigma_measure5_bar_star)+log((1/param_gamma_star)+(1/(1-param_gamma_star)))+log(1/param_sigmaa_star)+log(1/param_psiI_star)+log(1/param_sigmaL_star)+...
        log((1/param_aest_star)+(1/(1-param_aest_star)))+log((1/param_awest_star)+(1/(1-param_awest_star)))+log(1/param_BETA_EST_star);
    post_PF_star=prior_PF_star+lik_PF_star;
    r1_PF = exp(post_PF_star-post_PF+jac_PF-jac_PF_star);
    %r1_PF
    C1 = min(1,r1_PF);
    rng('shuffle')
    A_rand = rand;
    if A_rand<=C1
        A=A_star; 
        B=B_star;
        H=H_star;
        R=R_star;
        S2=S2_star;
        Phi=Phi_star;
        Phi11=Phi11_star;
        Phi12=Phi12_star;
        R00=R00_star;
        delta_ss=delta_ss_star;
        
        
        param_rhoR=param_rhoR_star;
        param_rhog=param_rhog_star;
        param_rhomu=param_rhomu_star;
        param_sigma_gest = param_sigma_gest_star;
        param_sigma_muest = param_sigma_muest_star;
        param_sigma_etaest = param_sigma_etaest_star;
        param_sigma_zest = param_sigma_zest_star;
        param_sigma_rest = param_sigma_rest_star;
        param_sigma_measure1_bar = param_sigma_measure1_bar_star;
        param_sigma_measure2_bar = param_sigma_measure2_bar_star;
        param_sigma_measure3_bar = param_sigma_measure3_bar_star;
        param_sigma_measure4_bar = param_sigma_measure4_bar_star;
        param_sigma_measure5_bar = param_sigma_measure5_bar_star;
        param_gammag = param_gammag_star;
        param_gammapi = param_gammapi_star;
        param_gammax = param_gammax_star;
        param_gamma = param_gamma_star;
        param_sigmaa = param_sigmaa_star;
        param_psip_bar = param_psip_bar_star;
        param_psiw_bar = param_psiw_bar_star;
        param_psiI = param_psiI_star;
        param_sigmaL = param_sigmaL_star;
        param_aest = param_aest_star;
        param_awest = param_awest_star;
        param_BETA_EST = param_BETA_EST_star;
        param_Gzest = param_Gzest_star;
        param_pibar_est = param_pibar_est_star;
        param_alpha = param_alpha_star;

        accept=accept+1;
        particles_eps=particles_eps_star;
        particles_first=particles_first_star;
        particles_second=particles_second_star;
        w=w_star;
        indx=indx_star;
        post_PF=post_PF_star;
        lik_indv=lik_indv_star;
        jac_PF=jac_PF_star;
        u_particles=u_particles_star;
        u_res=u_res_star;
    end
    
    end
    
    %end
    
    %else
       
    %keep the old ones    
        
        
        
        
    end
    
    
    for s=1:S
       [eps_ctraj(:,:,s),state_first_ctraj(:,:,s),state_second_ctraj(:,:,s)]=ancestor_trace_DSGE_2ndorder(particles_eps(:,:,:,s),particles_first(:,:,:,s),particles_second(:,:,:,s),w(:,:,s),indx(:,:,s));
    end
    
    
    theta=[logit_inverse(param_rhoR),logit_inverse(param_rhog),logit_inverse(param_rhomu),log(param_sigma_gest),...
        log(param_sigma_muest),log(param_sigma_etaest),log(param_sigma_zest),log(param_sigma_rest),...
        log(param_sigma_measure1_bar),log(param_sigma_measure2_bar),log(param_sigma_measure3_bar),log(param_sigma_measure4_bar),...
        log(param_sigma_measure5_bar),param_gammag,param_gammapi,param_gammax,logit_inverse(param_gamma),log(param_sigmaa),...
        param_psip_bar,param_psiw_bar,log(param_psiI),log(param_sigmaL),logit_inverse(param_aest),logit_inverse(param_awest),...
        log(param_BETA_EST),param_Gzest,param_pibar_est,param_alpha]; 

    thetasave(i,:)=theta;
 
    if i>=1000
       %scale=update_sigma(scale,C1,target_accept,i,D1);
       V1=cov(thetasave(1:i,:));
       %mean_param=mean(thetasave(1:i,:));
       V1=jitChol(V1);
    end
    %(i-adapt:end,1)
    if scale<=0.01
       scale=0.01;
     end
    
    Post.rhoR(i,1)=param_rhoR;
    Post.rhog(i,1)=param_rhog;
    Post.rhomu(i,1)=param_rhomu;
    Post.sigma_gest(i,1)=param_sigma_gest;
    Post.sigma_muest(i,1)=param_sigma_muest;
    Post.sigma_etaest(i,1)=param_sigma_etaest;
    Post.sigma_zest(i,1)=param_sigma_zest;
    Post.sigma_rest(i,1)=param_sigma_rest;
    Post.sigma_measure1(i,1)=param_sigma_measure1_bar/scale_param_shock;
    Post.sigma_measure2(i,1)=param_sigma_measure2_bar/scale_param_shock;
    Post.sigma_measure3(i,1)=param_sigma_measure3_bar/scale_param_shock;
    Post.sigma_measure4(i,1)=param_sigma_measure4_bar/scale_param_shock;
    Post.sigma_measure5(i,1)=param_sigma_measure5_bar/scale_param_shock;
    Post.gammag(i,1) = param_gammag;
    Post.gammapi(i,1) = param_gammapi;
    Post.gammax(i,1) = param_gammax;
    Post.gamma(i,1) = param_gamma;
    Post.sigmaa(i,1) = param_sigmaa;
    Post.psip(i,1) = param_psip_bar*scale_psip;
    Post.psiw(i,1) = param_psiw_bar*scale_psiw;
    Post.psiI(i,1) = param_psiI;
    Post.sigmaL(i,1) = param_sigmaL;
    Post.param_aest(i,1) = param_aest;
    Post.param_awest(i,1) = param_awest;
    Post.BETA_EST(i,1) = param_BETA_EST;
    Post.Gzest(i,1) = param_Gzest;
    Post.pibar_est(i,1) = param_pibar_est;
    Post.alpha(i,1) = param_alpha;
    
    
    
    %toc
    delete('PARAMSTORE.mat');
    if mod(i,250)==0
       save(['CorrPMMH','N',num2str(N),'S',num2str(S),'mediumScale_DataReal_disturbance.mat'],'Post');

    end
end
save(['CorrPMMH','N',num2str(N),'S',num2str(S),'mediumScale_DataReal_disturbance.mat'],'Post');

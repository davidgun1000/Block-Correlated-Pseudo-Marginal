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
%addpath('/home/davidg/dynare/lib/dynare/matlab')

%For KATANA
%addpath '/apps/dynare/4.4.3/lib/dynare/matlab'


%for GADI
addpath('/home/561/dg2271/dynare457/lib/dynare/matlab')
addpath('/home/561/dg2271/matio1.5.17') 
%parpool(10)
load('SmallScale_HS_data_full.mat');
yt=y_full(:,1:end);
param_zeta=0.2;
param_nu=1/11;
param_phi=200;




trim_first=0.25;
trim_second=0.75;

param_piA=4.01;
param_tau=0.90;
param_psi1=3.50;
param_psi2=1.12;
param_gammaQ=0.89;
param_rA=1.51;
param_BETA=1/(1+param_rA/400);

param_rhog=0.45;
param_rhoz=0.99;
param_rhoR=0.50;
param_sigma_r=0.04;
param_sigma_g= 0.09;
param_sigma_z=0.11;
param_sigma_measure1=0.60;
param_sigma_measure2=1.27;
param_sigma_measure3=0.04;


zeta=param_zeta;
nu=param_nu;
phi=param_phi;

piA=param_piA;
tau=param_tau;
psi1=param_psi1;
psi2=param_psi2;
gammaQ=param_gammaQ;
rA=param_rA;
BETA=param_BETA;

rhog=param_rhog;
rhoz=param_rhoz;
rhoR=param_rhoR;
sigma_r=param_sigma_r;
sigma_g=param_sigma_g;
sigma_z=param_sigma_z;
sigma_measure1=param_sigma_measure1;
sigma_measure2=param_sigma_measure2;
sigma_measure3=param_sigma_measure3;
%1-->gdp
eq_y=1;
%2--> inflation
eq_pi=2;
%3--> Interest Rate
eq_ffr=3;

save PARAMSTORE  BETA nu phi piA tau psi1 psi2 rhog rhoz rhoR sigma_z sigma_g sigma_r zeta gammaQ rA;

dynare eig_check noclearall nolog

B=B_true;
H=diag([sigma_measure1^2;sigma_measure2^2;sigma_measure3^2]);
dynare SmallScale_V4_2ndOrder_USE noclearall nolog 
f_11=[oo_.dr.ghx oo_.dr.ghu];
order=2;
optPruning=[];
[omega_ss,F1,F2,omega_ss_delta,F11,F22,F12,~,~,~,~,~,state_var_positions,control_var_postion] = statespaceform(optPruning,oo_,M_,f_11,order);
%clearvars -except sig order oo_ M_ omega_ss F1 F2 omega_ss_delta F11 F22 F12 F111 F112 F122 F222 F331 state_var_positions control_var_postion  BETTA DELTA ALFA RHO RHOE SIG sigma_a sigma_e...
%    H A B
prior_rhoR_mean=0.5;
prior_rhoR_std=0.2;
prior_rhog_mean=0.5;
prior_rhog_std=0.2;
prior_rhoz_mean=0.5;
prior_rhoz_std=0.2;
prior_sigma_r_mean=0.05;
prior_sigma_r_std=0.05;
prior_sigma_g_mean=0.05;
prior_sigma_g_std=0.05;
prior_sigma_z_mean=0.05;
prior_sigma_z_std=0.05;
prior_sigma_v0=10;
prior_sigma_s0=1;
prior_psi1_mean=2.41;
prior_psi1_std=0.5;
prior_psi2_mean=0.39;
prior_psi2_std=0.5;
prior_tau_mean=1.5;
prior_tau_std=0.5;
prior_piA_mean=3.1293;
prior_piA_std=0.5;
prior_gammaQ_mean=0.6075;
prior_gammaQ_std=0.5;
prior_rA_mean=0.7055;
prior_rA_std=0.5;
prior_sigma_measure_mean=2;
prior_sigma_measure_std=2;
prior_gamma_a=1;
prior_gamma_b=1;


Phi=F1;
R=F2;
S2=M_.Sigma_e;
Phi11=F11;
Phi12=F12;
R00=F22;
delta_ss=omega_ss_delta;

A(eq_y,1) = gammaQ;
A(eq_pi,1) = piA;
A(eq_ffr,1) = piA+rA+4*gammaQ;

S=100; %number of nodes
N=100;
burn=2000;
adapt=3000;
nit=25000;
nloop=burn+nit+adapt;
ne=size(S2,1);
[~, ns] = size(B);
T=size(yt,2);
num_prop=100;
N_prop=1000;
corr_param=0.9;
num_param=15;
optimal_scale=2.20^2/num_param;
num_mix=1;
corr_index=0;
disturbance_index=0;
w_prior=0.05;
w_rest=1-w_prior;

for s=1:num_prop
    [particles_eps_temp(:,:,:,s),particles_first_temp(:,:,:,s),particles_second_temp(:,:,:,s),w_temp(:,:,s),indx_temp(:,:,s),lik_temp(:,s)]=Mix_smc_SmallScale_2ndorder_bootstrap_alt(A,B,H,Phi,R,S2,Phi11,Phi12,R00,delta_ss,1000,yt,state_var_positions);
    [eps_ctraj_prop(:,:,s),~,~]=ancestor_trace_DSGE_2ndorder(particles_eps_temp(:,:,:,s),particles_first_temp(:,:,:,s),particles_second_temp(:,:,:,s),w_temp(:,:,s),indx_temp(:,:,s));
end

eps1_prop=(reshape(eps_ctraj_prop(1,:,:),T,num_prop))';
eps2_prop=(reshape(eps_ctraj_prop(2,:,:),T,num_prop))';
eps3_prop=(reshape(eps_ctraj_prop(3,:,:),T,num_prop))';

for j=1:T
    
    eps_collect=[eps1_prop(:,j),eps2_prop(:,j),eps3_prop(:,j)]; 
    covmat_CDKF(:,:,j)=cov(eps_collect);
    temp_var=cov(eps_collect);
    temp_mean=mean(eps_collect);
    temp_w=1;


    temp_w=temp_w.*w_rest;
    w_mix(:,j)=[w_prior;temp_w];
    mean_mix(:,:,j)=[zeros(1,ne);temp_mean];
    covmat_mix(:,:,:,j)=cat(3,eye(ne),temp_var(:,:,1:num_mix));
    
end 

%-------------------
%generating initial proposal

for s=1:num_prop
   for t=1:T
       unif_mix_prop(1,1:N_prop)=sort(rand(1,N_prop));
       for k=1:num_mix+1
           id_prop=(unif_mix_prop(1,1:N_prop)>=sum(w_mix(1:(k-1),t))) & (unif_mix_prop(1,1:N_prop)<=sum(w_mix(1:k,t)));
           id_mix_prop(t,id_prop,s)=k;
           
       end
   end
   u_particles_prop(:,1:N_prop,:,s)=randn(ne,N_prop,T);
   u_res_prop(:,1:N_prop,s)=rand(T,N_prop);
end
A_prop=A;
B_prop=B;
H_prop=H;
Phi_prop=F1;
R_prop=F2;
S2_prop=M_.Sigma_e;
Phi11_prop=F11;
Phi12_prop=F12;
R00_prop=F22;
delta_ss_prop=omega_ss_delta; 
parfor s=1:num_prop
   [particles_eps_prop(:,:,:,s),particles_first_prop(:,:,:,s),particles_second_prop(:,:,:,s),w_prop(:,:,s),indx_prop(:,:,s),lik_prop(:,s)]=Mix_smc_SmallScale_2ndorder_PMMH_mix(A_prop,B_prop,H_prop,Phi_prop,R_prop,S2_prop,Phi11_prop,Phi12_prop,...
        R00_prop,delta_ss_prop,N_prop,yt,u_particles_prop(:,:,:,s),...
        u_res_prop(:,:,s),...
        w_mix,mean_mix,covmat_mix,...
        id_mix_prop(:,:,s),state_var_positions,0);
   [eps_ctraj_prop(:,:,s),~,~]=ancestor_trace_DSGE_2ndorder(particles_eps_prop(:,:,:,s),particles_first_prop(:,:,:,s),particles_second_prop(:,:,:,s),w_prop(:,:,s),indx_prop(:,:,s));
end

eps1_prop=(reshape(eps_ctraj_prop(1,:,:),T,num_prop))';
eps2_prop=(reshape(eps_ctraj_prop(2,:,:),T,num_prop))';
eps3_prop=(reshape(eps_ctraj_prop(3,:,:),T,num_prop))';

for j=1:T
    
    eps_collect=[eps1_prop(:,j),eps2_prop(:,j),eps3_prop(:,j)]; 
    covmat_CDKF(:,:,j)=cov(eps_collect);
    temp_var=cov(eps_collect);
    temp_mean=mean(eps_collect);
    temp_w=1;
    
     temp_w=temp_w.*w_rest;
     w_mix(:,j)=[w_prior;temp_w];
     mean_mix(:,:,j)=[zeros(1,ne);temp_mean];
     covmat_mix(:,:,:,j)=cat(3,eye(ne),temp_var(:,:,1:num_mix));


end 


for s=1:S
    for t=1:T
       unif_mix(1,1:N)=sort(rand(1,N));
       for k=1:num_mix+1
           id=(unif_mix(1,1:N)>=sum(w_mix(1:(k-1),t))) & (unif_mix(1,1:N)<=sum(w_mix(1:k,t)));
           id_mix(t,id,s)=k;
           
       end
    end
    
   
    
    u_particles(:,1:N,:,s)=randn(ne,N,T);
    u_res(:,1:N,s)=rand(T,N);
    
    
    
end
parfor s=1:S
    [particles_eps(:,:,:,s),particles_first(:,:,:,s),particles_second(:,:,:,s),w(:,:,s),indx(:,:,s),lik_indv(:,s)]=Mix_smc_SmallScale_2ndorder_PMMH_mix(A,B,H,Phi,R,S2,Phi11,Phi12,R00,delta_ss,N,yt,u_particles(:,:,:,s),u_res(:,:,s),...
        w_mix,mean_mix,covmat_mix,...
        id_mix(:,:,s),state_var_positions,corr_index,disturbance_index);
    [eps_ctraj(:,:,s),state_first_ctraj(:,:,s),state_second_ctraj(:,:,s)]=ancestor_trace_DSGE_2ndorder(particles_eps(:,:,:,s),particles_first(:,:,:,s),particles_second(:,:,:,s),w(:,:,s),indx(:,:,s));
end
lik_indv_appr1 = sum(lik_indv);
lik_indv_appr1_sort = sort(lik_indv_appr1);
lik_indv_appr1_trim = lik_indv_appr1_sort((round(trim_first*S):round(trim_second*S)));
maxw_appr1 = max(lik_indv_appr1_trim);
lik_PF = log(mean(exp(lik_indv_appr1_trim-maxw_appr1)))+maxw_appr1; % average of the log likelihood

prior_PF=log_truncatednormal(param_tau,prior_tau_mean,prior_tau_std,0,inf)+...
      log_truncatednormal(param_psi1,prior_psi1_mean,prior_psi1_std,0,inf)+log_truncatednormal(param_psi2,prior_psi2_mean,prior_psi2_std,0,inf)+...
      log_truncatednormal(param_rhoR,prior_rhoR_mean,prior_rhoR_std,0,1)+...
      log_truncatednormal(param_rhog,prior_rhog_mean,prior_rhog_std,0,1)+log_truncatednormal(param_rhoz,prior_rhoz_mean,prior_rhoz_std,0,1)+...
      log_IG_PDF_used(param_sigma_r,prior_sigma_v0/2,prior_sigma_s0/2)+log_IG_PDF_used(param_sigma_g,prior_sigma_v0/2,prior_sigma_s0/2)+...
      log_IG_PDF_used(param_sigma_z,prior_sigma_v0/2,prior_sigma_s0/2)+...
      log_truncatednormal(param_rA,prior_rA_mean,prior_rA_std,0,inf)+...
      log_truncatednormal(param_piA,prior_piA_mean,prior_piA_std,0,inf)+log(normpdf(param_gammaQ,prior_gammaQ_mean,prior_gammaQ_std))+...
      log(gampdf(param_sigma_measure1,prior_gamma_a,prior_gamma_b))+...
      log(gampdf(param_sigma_measure2,prior_gamma_a,prior_gamma_b))+...
      log(gampdf(param_sigma_measure3,prior_gamma_a,prior_gamma_b));
      
post_PF=prior_PF+lik_PF;

jac_PF=log(1/param_tau)+log(1/param_psi1)+log(1/param_psi2)+...
    log((1/param_rhoR)+(1/(1-param_rhoR)))+log((1/param_rhog)+(1/(1-param_rhog)))+log((1/param_rhoz)+(1/(1-param_rhoz)))+...
    log(1/param_sigma_r)+log(1/param_sigma_g)+log(1/param_sigma_z)+log(1/param_piA)+log(1/param_rA)+...
    log(1/param_sigma_measure1)+log(1/param_sigma_measure2)+log(1/param_sigma_measure3);

load('covmat_smallscale_2ndorder_nonfixed.mat');
D1=15;
V1=0.01*eye(D1);
scale=0.01;
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
    rng('shuffle')
    theta=[log(param_tau),log(param_psi1),log(param_psi2),logit_inverse(param_rhoR),...
           logit_inverse(param_rhog),logit_inverse(param_rhoz),log(param_sigma_r),log(param_sigma_g),log(param_sigma_z),...
           log(param_rA),log(param_piA),param_gammaQ,log(param_sigma_measure1),log(param_sigma_measure2),log(param_sigma_measure3)]; 
            
    rng('shuffle')
    R1=mvnrnd(theta,scale.*V1);
    
    param_tau_star=exp(R1(1,1));
    param_psi1_star=exp(R1(1,2));
    param_psi2_star=exp(R1(1,3));
    param_rhoR_star=logitcdf(R1(1,4));
    param_rhog_star=logitcdf(R1(1,5));
    param_rhoz_star=logitcdf(R1(1,6));
    param_sigma_r_star=exp(R1(1,7));
    param_sigma_g_star=exp(R1(1,8));
    param_sigma_z_star=exp(R1(1,9));
    param_rA_star=exp(R1(1,10));
    param_piA_star=exp(R1(1,11));
    param_gammaQ_star=R1(1,12);
    param_sigma_measure1_star=exp(R1(1,13));
    param_sigma_measure2_star=exp(R1(1,14));
    param_sigma_measure3_star=exp(R1(1,15));
    
    param_BETA_star=1/(1+param_rA_star/400);
    
    zeta=param_zeta;
    nu=param_nu;
    phi=param_phi;

    piA=param_piA_star;
    tau=param_tau_star;
    psi1=param_psi1_star;
    psi2=param_psi2_star;
    gammaQ=param_gammaQ_star;
    rA=param_rA_star;
    BETA=param_BETA_star;
    rhog=param_rhog_star;
    rhoz=param_rhoz_star;
    rhoR=param_rhoR_star;
    sigma_r=param_sigma_r_star;
    sigma_g=param_sigma_g_star;
    sigma_z=param_sigma_z_star;
    sigma_measure1=param_sigma_measure1_star;
    sigma_measure2=param_sigma_measure2_star;
    sigma_measure3=param_sigma_measure3_star;

    save PARAMSTORE  BETA nu phi piA tau psi1 psi2 rhog rhoz rhoR sigma_z sigma_g sigma_r zeta gammaQ rA;
    
    eig_check
    
    if size(find(real(abs(oo_.dr.eigval))>1),1)==M_.nsfwrd
    
    
    SmallScale_V4_2ndOrder_USE
        
    f_11=[oo_.dr.ghx oo_.dr.ghu];
    order=2;
    optPruning=[];
    [omega_ss,F1,F2,omega_ss_delta,F11,F22,F12,~,~,~,~,~,state_var_positions,control_var_postion] = statespaceform(optPruning,oo_,M_,f_11,order);

    A_star(eq_y,1) = gammaQ;
    A_star(eq_pi,1) = piA;
    A_star(eq_ffr,1) = piA+rA+4*gammaQ;
    %A_star=A;
    B_star=B;
    H_star=diag([sigma_measure1^2;sigma_measure2^2;sigma_measure3^2]);
    Phi_star=F1;
    R_star=F2;
    S2_star=M_.Sigma_e;
    Phi11_star=F11;
    Phi12_star=F12;
    R00_star=F22;
    delta_ss_star=omega_ss_delta;
    
    rng('shuffle')
    
    for s=1:S
        for t=1:T
           unif_mix(1,1:N)=sort(rand(1,N));
           for k=1:num_mix+1
               id=(unif_mix(1,1:N)>=sum(w_mix(1:(k-1),t))) & (unif_mix(1,1:N)<=sum(w_mix(1:k,t)));
               id_mix_star(t,id,s)=k;

           end
        end
    
   
    
    u_particles_star(:,1:N,:,s)=randn(ne,N,T);
    u_res_star(:,1:N,s)=rand(T,N);
    
    
    
    end
               
    %tic
    parfor s=1:S
        [particles_eps_star(:,:,:,s),particles_first_star(:,:,:,s),particles_second_star(:,:,:,s),w_star(:,:,s),indx_star(:,:,s),lik_indv_star(:,s)]=Mix_smc_SmallScale_2ndorder_PMMH_mix(A_star,B_star,H_star,Phi_star,R_star,S2_star,...
            Phi11_star,Phi12_star,R00_star,delta_ss_star,N,yt,u_particles_star(:,:,:,s),u_res_star(:,:,s),...
            w_mix,mean_mix,covmat_mix,...
            id_mix_star(:,:,s),state_var_positions,corr_index,disturbance_index);        
    end
    %toc
    %time(i,1)=toc;
    if sum(isnan(lik_indv_star))>0 | sum(sum(sum(sum(isnan(particles_eps_star),4),3),2))>0 
       A=A; 
       B=B;
       H=H;
       R=R;
       S2=S2;
       Phi=Phi;
       Phi11=Phi11;
       Phi12=Phi12;
       R00=R00;
       delta_ss=delta_ss;
       param_piA=param_piA;
       param_tau=param_tau;
       param_psi1=param_psi1;
       param_psi2=param_psi2;
       param_gammaQ=param_gammaQ;
       param_rA=param_rA;
       param_BETA=param_BETA;
       param_rhog=param_rhog;
       param_rhoz=param_rhoz;
       param_rhoR=param_rhoR;
       param_sigma_r=param_sigma_r;
       param_sigma_g=param_sigma_g;
       param_sigma_z=param_sigma_z;
       param_sigma_measure1=param_sigma_measure1;
       param_sigma_measure2=param_sigma_measure2;
       param_sigma_measure3=param_sigma_measure3;

       particles_eps=particles_eps;
       particles_first=particles_first;
       particles_second=particles_second;
       w=w;
       indx=indx;
       post_PF=post_PF;
       %post_CDKF=post_CDKF;
       %jac_CDKF=jac_CDKF;
       jac_PF=jac_PF;
       
       u_particles=u_particles;
       u_res=u_res;
       id_mix=id_mix; 
       lik_indv=lik_indv;
       
       
       
    else
    
    lik_indv_appr1_star = sum(lik_indv_star);
    lik_indv_appr1_sort_star = sort(lik_indv_appr1_star);
    lik_indv_appr1_trim_star = lik_indv_appr1_sort_star((round(trim_first*S):round(trim_second*S)));
    maxw_appr1_star = max(lik_indv_appr1_trim_star);
    lik_PF_star = log(mean(exp(lik_indv_appr1_trim_star-maxw_appr1_star)))+maxw_appr1_star; % average of the log likelihood
    %interacting step
    
    prior_PF_star=log_truncatednormal(param_tau_star,prior_tau_mean,prior_tau_std,0,inf)+...
      log_truncatednormal(param_psi1_star,prior_psi1_mean,prior_psi1_std,0,inf)+log_truncatednormal(param_psi2_star,prior_psi2_mean,prior_psi2_std,0,inf)+...
      log_truncatednormal(param_rhoR_star,prior_rhoR_mean,prior_rhoR_std,0,1)+...
      log_truncatednormal(param_rhog_star,prior_rhog_mean,prior_rhog_std,0,1)+log_truncatednormal(param_rhoz_star,prior_rhoz_mean,prior_rhoz_std,0,1)+...
      log_IG_PDF_used(param_sigma_r_star,prior_sigma_v0/2,prior_sigma_s0/2)+log_IG_PDF_used(param_sigma_g_star,prior_sigma_v0/2,prior_sigma_s0/2)+...
      log_IG_PDF_used(param_sigma_z_star,prior_sigma_v0/2,prior_sigma_s0/2)+...
      log_truncatednormal(param_rA_star,prior_rA_mean,prior_rA_std,0,inf)+...
      log_truncatednormal(param_piA_star,prior_piA_mean,prior_piA_std,0,inf)+log(normpdf(param_gammaQ_star,prior_gammaQ_mean,prior_gammaQ_std))+...
      log(gampdf(param_sigma_measure1_star,prior_gamma_a,prior_gamma_b))+...
      log(gampdf(param_sigma_measure2_star,prior_gamma_a,prior_gamma_b))+...
      log(gampdf(param_sigma_measure3_star,prior_gamma_a,prior_gamma_b));
    
    jac_PF_star=log(1/param_tau_star)+log(1/param_psi1_star)+log(1/param_psi2_star)+...
        log((1/param_rhoR_star)+(1/(1-param_rhoR_star)))+log((1/param_rhog_star)+(1/(1-param_rhog_star)))+log((1/param_rhoz_star)+(1/(1-param_rhoz_star)))+...
        log(1/param_sigma_r_star)+log(1/param_sigma_g_star)+log(1/param_sigma_z_star)+log(1/param_piA_star)+log(1/param_rA_star)+...
        log(1/param_sigma_measure1_star)+log(1/param_sigma_measure2_star)+log(1/param_sigma_measure3_star);
    

      
    post_PF_star=prior_PF_star+lik_PF_star;
    r1_PF = exp(post_PF_star-post_PF+jac_PF-jac_PF_star);
    
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
        param_piA=param_piA_star;
        param_tau=param_tau_star;
        param_psi1=param_psi1_star;
        param_psi2=param_psi2_star;
        param_gammaQ=param_gammaQ_star;
        param_rA=param_rA_star;
        param_BETA=param_BETA_star;
        param_rhog=param_rhog_star;
        param_rhoz=param_rhoz_star;
        param_rhoR=param_rhoR_star;
        param_sigma_r=param_sigma_r_star;
        param_sigma_g=param_sigma_g_star;
        param_sigma_z=param_sigma_z_star;
        param_sigma_measure1=param_sigma_measure1_star;
        param_sigma_measure2=param_sigma_measure2_star;
        param_sigma_measure3=param_sigma_measure3_star;

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
        id_mix=id_mix_star;
    end
    
    end
    
    %end
    
    else
        A=A; 
       B=B;
       H=H;
       R=R;
       S2=S2;
       Phi=Phi;
       Phi11=Phi11;
       Phi12=Phi12;
       R00=R00;
       delta_ss=delta_ss;
       param_piA=param_piA;
       param_tau=param_tau;
       param_psi1=param_psi1;
       param_psi2=param_psi2;
       param_gammaQ=param_gammaQ;
       param_rA=param_rA;
       param_BETA=param_BETA;
       param_rhog=param_rhog;
       param_rhoz=param_rhoz;
       param_rhoR=param_rhoR;
       param_sigma_r=param_sigma_r;
       param_sigma_g=param_sigma_g;
       param_sigma_z=param_sigma_z;
       param_sigma_measure1=param_sigma_measure1;
        param_sigma_measure2=param_sigma_measure2;
        param_sigma_measure3=param_sigma_measure3;

       particles_eps=particles_eps;
       particles_first=particles_first;
       particles_second=particles_second;
       w=w;
       indx=indx;
       post_PF=post_PF;
       %post_CDKF=post_CDKF;
       %jac_CDKF=jac_CDKF;
       jac_PF=jac_PF;
       
       u_particles=u_particles;
       u_res=u_res;
       id_mix=id_mix; 
       lik_indv=lik_indv;
        
        
        
        
        
        
    end
    
    
    for s=1:S
       [eps_ctraj(:,:,s),state_first_ctraj(:,:,s),state_second_ctraj(:,:,s)]=ancestor_trace_DSGE_2ndorder(particles_eps(:,:,:,s),particles_first(:,:,:,s),particles_second(:,:,:,s),w(:,:,s),indx(:,:,s));
    end
    
    
    theta=[log(param_tau),log(param_psi1),log(param_psi2),logit_inverse(param_rhoR),...
            logit_inverse(param_rhog),logit_inverse(param_rhoz),log(param_sigma_r),log(param_sigma_g),log(param_sigma_z),...
            log(param_rA),log(param_piA),gammaQ,log(param_sigma_measure1),log(param_sigma_measure2),log(param_sigma_measure3)];

    thetasave(i,:)=theta;
 
    if i>=burn
       scale=update_sigma(scale,C1,target_accept,i,D1);
       V1=cov(thetasave(1:i,:));
       %mean_param=mean(thetasave(1:i,:));
       V1=jitChol(V1);
    end
    %(i-adapt:end,1)
    if scale<=0.01
       scale=0.01;
     end
    eps1_store=(reshape(eps_ctraj(1,:,:),T,S))';
    eps2_store=(reshape(eps_ctraj(2,:,:),T,S))';
    eps3_store=(reshape(eps_ctraj(3,:,:),T,S))';
    
    
    for j=1:T
    
        eps_collect=[eps1_store(:,j),eps2_store(:,j),eps3_store(:,j)]; 
        covmat_CDKF(:,:,j)=cov(eps_collect);
        temp_var=cov(eps_collect);
        temp_mean=mean(eps_collect);
        temp_w=1;


        temp_w=temp_w.*w_rest;
        w_mix(:,j)=[w_prior;temp_w];
        mean_mix(:,:,j)=[zeros(1,ne);temp_mean];
        covmat_mix(:,:,:,j)=cat(3,eye(ne),temp_var(:,:,1:num_mix));
    
    end     
    
    Post.piA(i,1)=param_piA;
    Post.tau(i,1)=param_tau;
    Post.psi1(i,1)=param_psi1;
    Post.psi2(i,1)=param_psi2;
    Post.gammaQ(i,1)=param_gammaQ;
    Post.rA(i,1)=param_rA;
    Post.rhoR(i,1)=param_rhoR;
    Post.rhog(i,1)=param_rhog;
    Post.rhoz(i,1)=param_rhoz;
    Post.sigma_r(i,1)=param_sigma_r;
    Post.sigma_g(i,1)=param_sigma_g;
    Post.sigma_z(i,1)=param_sigma_z;
    Post.var(i,1)=var(sum(lik_indv));
    Post.sigma_measure1(i,1)=param_sigma_measure1;
    Post.sigma_measure2(i,1)=param_sigma_measure2;
    Post.sigma_measure3(i,1)=param_sigma_measure3;

    %toc
    delete('PARAMSTORE.mat');
    if mod(i,250)==0
       save(['ParPMMH_appr1_nocorr_',num2str(trim_first),num2str(trim_second),'N',num2str(N),'S',num2str(S),'SmallScale_Data1_disturbance.mat'],'Post');

    end
end
save(['ParPMMH_appr1_nocorr_',num2str(trim_first),num2str(trim_second),'N',num2str(N),'S',num2str(S),'SmallScale_Data1_disturbance.mat'],'Post');

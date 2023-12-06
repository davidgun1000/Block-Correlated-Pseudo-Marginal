

clear
clc
close all
delete *.asv
tic
%addpath 'C:\dynare\4.6.1\matlab'
%addpath 'C:\Users\dgunawan\Dropbox\Gunawan_Kohn\DavidCodev1\RBC_sim'
%addpath 'C:\dynare\4.4.3\matlab'
%addpath('C:\dynare\4.5.7\matlab');
%addpath 'D:\Dropbox\Gunawan_Kohn\DavidCodev1\RBC_sim'
addpath '/Applications/Dynare/4.5.7/matlab'
%addpath('/home/davidg/dynare/lib/dynare/matlab')

load('sim_SMALL_SCALE_T80.mat');
yt=y_simul;

param_zeta=zeta_true;
param_nu=nu_true;
param_BETA=BETA_true;
param_phi=phi_true;
param_piss=piss_true;
param_gamma=gamma_true;

param_tau=tau_true;
param_psi1=psi1_true;
param_psi2=psi2_true;
param_rhog=rhog_true;
param_rhoz=rhoz_true;
param_rhoR=rhoR_true;
param_sigma_r=sigma_r_true;
param_sigma_g=sigma_g_true;
param_sigma_z=sigma_z_true;

zeta=param_zeta;
nu=param_nu;
BETA=param_BETA;
phi=param_phi;
piss=param_piss;
gamma=param_gamma;

tau=param_tau;
psi1=param_psi1;
psi2=param_psi2;
rhog=param_rhog;
rhoz=param_rhoz;
rhoR=param_rhoR;
sigma_r=param_sigma_r;
sigma_g=param_sigma_g;
sigma_z=param_sigma_z;

save PARAMSTORE BETA nu phi piss tau psi1 psi2 rhog rhoz rhoR sigma_z sigma_g sigma_r zeta gamma;

B=B_true;
H=H_true;
dynare SmallScale noclearall nolog 
f_11=[oo_.dr.ghx oo_.dr.ghu];
order=2;
optPruning=[];
[omega_ss,F1,F2,omega_ss_delta,F11,F22,F12,~,~,~,~,~,state_var_positions,control_var_postion] = statespaceform(optPruning,oo_,M_,f_11,order);

Phi=F1;
R=F2;
S2=M_.Sigma_e;
Phi11=F11;
Phi12=F12;
R00=F22;
delta_ss=omega_ss_delta;

S=100; %number of nodes
group_size=[1];
group_num=length(group_size);
num_block=S./group_size;
N=250;
burn=5000;
nit=20000;
nloop=burn+nit;
ne=size(S2,1);
[~, ns] = size(B);
T=size(yt,2);
%load('init_prop.mat');
%parpool(12)
%parpool(28)
num_prop=100;
N_prop=2000;
corr_param=0.9999;
num_param=4;
optimal_scale=2.16^2/num_param;
num_mix=1;
w_prior=0.05;
w_rest=1-w_prior;

parfor s=1:num_prop
    [particles_eps_temp(:,:,:,s),particles_first_temp(:,:,:,s),particles_second_temp(:,:,:,s),w_temp(:,:,s),indx_temp(:,:,s),lik_temp(s,1)]=Mix_smc_RBC_2ndorder_bootstrap_alt(B,H,Phi,R,S2,Phi11,Phi12,R00,delta_ss,N_prop,yt,state_var_positions);
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
   [particles_eps_prop(:,:,:,s),particles_first_prop(:,:,:,s),particles_second_prop(:,:,:,s),w_prop(:,:,s),indx_prop(:,:,s),lik_prop(s,1)]=Mix_smc_RBC_2ndorder_PMMH_mix(B_prop,H_prop,Phi_prop,R_prop,S2_prop,Phi11_prop,Phi12_prop,...
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

save('prop_T80_SmallScale.mat','mean_mix','covmat_mix');
%checking correlation
load('simdata5_100.mat');
dim_y=size(yt,1);
D=0;
C=0;
M=eye(dim_y);
fs_init = zeros(dim_y,1);
omega_fs_init = eye(dim_y);
G=G_true;

S=[1000];
N=[100];

group_size = 1;
num_block=S./group_size;
corr_param=0.9999;

length_S=length(S);
reps=100;
for d=1:length_S

for r=1:reps
    r
for s=1:S(d,1)
    u_particles_init(:,:,s) = randn(dim_y,N);
    u_particles(:,:,:,s)=randn(dim_y,N,T);
    u_res(:,:,s)=rand(T,N);
end        
        
parfor s=1:S
    [~,~,~,lik_indv(:,s)]=bootstrap_LGSS_corr(yt,D,Z,sigma_u,C,G,M,sigma_e,u_particles_init(:,:,s),u_particles(:,:,:,s),u_res(:,:,s),N);
end
% temp1_appr3 = mean(lik_indv,2);
% temp2_appr3 = sum(lik_indv);
% sig2_appr3 = var(temp2_appr3');
% lik_appr3(r,1) = sum(temp1_appr3)+((sig2_appr3/2)*(1-(1/S)));




for t=1:T
    maxw_appr2 = max(lik_indv(t,:));
    lik_indv_appr2(t,:) = log(mean(exp(lik_indv(t,:)-maxw_appr2)))+maxw_appr2;
end
lik_appr2(r,1) = sum(lik_indv_appr2);



u_particles_init_star = u_particles_init;
u_particles_star=u_particles;
u_res_star=u_res;

k_group=randperm(S(d,1),1);
group_ind=(k_group-1)*group_size+1:k_group*group_size;

for s=group_ind(1):group_ind(end)
    u_particles_init_star(:,:,s) = corr_param*u_particles_init(:,:,s)+sqrt(1-corr_param^2)*randn(dim_y,N);    
    u_particles_star(:,:,:,s)=corr_param*u_particles(:,:,:,s)+sqrt(1-corr_param^2)*randn(dim_y,N,T);
    u_res_norm=norminv(u_res(:,:,s));
    u_res_norm_star=corr_param*u_res_norm+sqrt(1-corr_param^2)*randn(T,N);
    u_res_star(:,:,s)=normcdf(u_res_norm_star);
end    
    
parfor s=1:S(d,1)
    [~,~,~,lik_indv_star(:,s)]=bootstrap_LGSS_corr(yt,D,Z,sigma_u,C,G,M,sigma_e,u_particles_init(:,:,s),u_particles_star(:,:,:,s),u_res_star(:,:,s),N);
end
% temp1_appr3_star = mean(lik_indv_star,2);
% temp2_appr3_star = sum(lik_indv_star);
% sig2_appr3_star = var(temp2_appr3_star');
% lik_appr3_star(r,1) = sum(temp1_appr3_star)+((sig2_appr3_star/2)*(1-(1/S)));



for t=1:T
    maxw_appr2_star = max(lik_indv_star(t,:));
    lik_indv_appr2_star(t,:) = log(mean(exp(lik_indv_star(t,:)-maxw_appr2_star)))+maxw_appr2_star;
end
lik_appr2_star(r,1) = sum(lik_indv_appr2_star);    
    
    
    
    
    
    
end
corr_loglik(d,1)=corr(lik_appr2,lik_appr2_star);


end
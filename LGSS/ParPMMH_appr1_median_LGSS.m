
parpool(10)
load('simdata10_300.mat');
T = size(yt,2);
burn=5000;
nit=20000;
nloop=burn+nit;
dim_y=size(yt,1);
D=0;
C=0;
M=eye(dim_y);
G=G_true;
S=101;
N=250;

psi = 0.4;
for k=1:dim_y
    for j=1:k
        G(k,j)=psi^(abs(k-j)+1);
    end
end
prior_psi = log(unifpdf(psi,-1,1));
for s=1:S
    u_particles_init(:,:,s) = randn(dim_y,N);
    u_particles(:,:,:,s)=randn(dim_y,N,T);
    u_res(:,:,s)=rand(T,N);
end

parfor s=1:S
    [~,~,~,lik_indv(:,s)]=bootstrap_LGSS_corr(yt,D,Z,sigma_u,C,G,M,sigma_e,u_particles_init(:,:,s),u_particles(:,:,:,s),u_res(:,:,s),N);
end

lik_indv_appr1 = sum(lik_indv);
lik_indv_appr1_sort = sort(lik_indv_appr1);
lik_indv_appr1_median = lik_indv_appr1_sort((S+1)/2);
lik_appr1 = lik_indv_appr1_median;







post = lik_appr1 + prior_psi;
jac_PF=log((1/psi)+(1/(1-psi)));   
accept=0;
D1=1;
scale = 0.01;
target_accept=0.2;

V1=eye(D1);

group_size = 1;
num_block=S./group_size;
corr_param=0.9;





tic
for i=1:nloop
    %i
    %accept
    %psi
    A1 = rand();
    theta = logit_inverse(psi);
    R1 = mvnrnd(theta,scale*V1);
    psi_star = logit_cdf_com(R1(1,1));
    for k=1:dim_y
        for j=1:k
            G_star(k,j)=psi_star^(abs(k-j)+1);
        end
    end
    prior_psi_star = log(unifpdf(psi_star,-1,1));
    
    u_particles_init_star = u_particles_init;
    u_particles_star=u_particles;
    u_res_star=u_res;

    k_group=randperm(S,1);
    group_ind=(k_group-1)*group_size+1:k_group*group_size;

    for s=group_ind(1):group_ind(end)
        u_particles_init_star(:,:,s) = corr_param*u_particles_init(:,:,s)+sqrt(1-corr_param^2)*randn(dim_y,N);    
        u_particles_star(:,:,:,s)=corr_param*u_particles(:,:,:,s)+sqrt(1-corr_param^2)*randn(dim_y,N,T);
        u_res_norm=norminv(u_res(:,:,s));
        u_res_norm_star=corr_param*u_res_norm+sqrt(1-corr_param^2)*randn(T,N);
        u_res_star(:,:,s)=normcdf(u_res_norm_star);
    end    
    parfor s=1:S
        [~,~,~,lik_indv_star(:,s)]=bootstrap_LGSS_corr(yt,D,Z,sigma_u,C,G_star,M,sigma_e,u_particles_init_star(:,:,s),u_particles_star(:,:,:,s),u_res_star(:,:,s),N);
    end
    lik_indv_appr1_star = sum(lik_indv_star);
    lik_indv_appr1_sort_star = sort(lik_indv_appr1_star);
    lik_indv_appr1_median_star = lik_indv_appr1_sort_star((S+1)/2);
    lik_appr1_star = lik_indv_appr1_median_star;
    
    
    post_star = lik_appr1_star+prior_psi_star;
    jac_PF_star=log((1/psi_star)+(1/(1-psi_star)));   
    r1_PF = exp(post_star-post+jac_PF-jac_PF_star);
    C1 = min(1,r1_PF);
    if A1<=C1
       G=G_star;
       psi=psi_star;
       post=post_star;
       accept=accept+1; 
       jac_PF=jac_PF_star;  
       u_particles_init=u_particles_init_star;
       u_particles=u_particles_star;
       u_res=u_res_star;
    end
    thetasave(i,:)=[logit_inverse(psi)];
 
%     if i>=burn
%         V1=cov(thetasave(1:i,:));
%         V1=jitChol(V1);
%     end
    scale=update_sigma(scale,C1,target_accept,i,1);
    Post.psi(i,1) = psi;
    if mod(i,1000)==0
       save(['ParPMMH_appr1_median','dim',num2str(dim_y),'T',num2str(T),'N',num2str(N),'S',num2str(S),'.mat'],'Post');
 
    end
    
end
cputime=toc;
save(['ParPMMH_appr1_median','dim',num2str(dim_y),'T',num2str(T),'N',num2str(N),'S',num2str(S),'.mat'],'Post');




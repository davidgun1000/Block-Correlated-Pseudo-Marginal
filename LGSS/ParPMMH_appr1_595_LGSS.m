
%parpool(8)
load('simdata10_200.mat');
T = size(yt,2);
burn=5000;
nit=20000;
nloop=burn+nit;
dim_y=size(yt,1);
D=0;
C=0;
M=eye(dim_y);
G=G_true;
S=100;
N=100;

trim_first = 0.05;
trim_second = 0.95;

psi = 0.4;
for k=1:dim_y
    for j=1:k
        G(k,j)=psi^(abs(k-j)+1);
    end
end
prior_psi = log(unifpdf(psi,-1,1));
for s=1:S
    s
    u_particles_init(:,:,s) = randn(dim_y,N);
    u_particles(:,:,:,s)=randn(dim_y,N,T);
    u_res(:,:,s)=rand(T,N);
end

parfor s=1:S
    s
    [~,~,~,lik_indv(:,s)]=bootstrap_LGSS_corr(yt,D,Z,sigma_u,C,G,M,sigma_e,u_particles_init(:,:,s),u_particles(:,:,:,s),u_res(:,:,s),N);
end

lik_indv_appr1 = sum(lik_indv);
lik_indv_appr1_sort = sort(lik_indv_appr1);
lik_indv_appr1_trim = lik_indv_appr1_sort((round(trim_first*S):round(trim_second*S)));
maxw_appr1 = max(lik_indv_appr1_trim);
lik_appr1 = log(mean(exp(lik_indv_appr1_trim-maxw_appr1)))+maxw_appr1; % average of the log likelihood
post = lik_appr1 + prior_psi;
jac_PF=log((1/psi)+(1/(1-psi)));   
accept=0;
D1=1;
V1=eye(D1);
scale=0.01;
target_accept=0.2;
group_size = 1;
num_block=S./group_size;
corr_param=0.9;





for i=1:nloop
    tic

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
    lik_indv_appr1_trim_star = lik_indv_appr1_sort_star((round(trim_first*S):round(trim_second*S)));
    maxw_appr1_star = max(lik_indv_appr1_trim_star);
    lik_appr1_star = log(mean(exp(lik_indv_appr1_trim_star-maxw_appr1_star)))+maxw_appr1_star; % average of the log likelihood
    
    
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
 
    %if i>=burn
    %    V1=cov(thetasave(1:i,:));
    %    V1=jitChol(V1);
    %end
    scale=update_sigma(scale,C1,target_accept,i,1);
    Post.psi(i,1) = psi;
    if mod(i,1000)==0
       save(['ParPMMH_appr1_',num2str(trim_first),num2str(trim_second),'dim',num2str(dim_y),'T',num2str(T),'N',num2str(N),'S',num2str(S),'.mat'],'Post');
 
    end
    toc
end
cputime=toc;
save(['ParPMMH_appr1_',num2str(trim_first),num2str(trim_second),'dim',num2str(dim_y),'T',num2str(T),'N',num2str(N),'S',num2str(S),'.mat'],'Post','cputime');



% A=A_true;
% vech_A=vech(A);
% corr_param=0.99;
% num_param = length(vech_A);
% num_mix=1;
% 
% corr_index=1;
% disturbance_index=1;
% [~, ns] = size(A);
% ne=ns;
% w_prior = 0.05;
% w_rest = 1 - w_prior;
% 
% %parpool(28)
% num_prop=100;
% 
%  
% for s=1:S
%      u_particles(:,1:N,:,s)=randn(ne,N,T);
%      u_res(:,1:N,s)=rand(T,N);
% end
%  
% for s=1:S    
%     [particles_eps_temp(:,:,:,s),w_temp(:,:,s),indx_temp(:,:,s),lik_temp(s,1)]=bootstrap_LGSS_PMMH(B,D,A,C,5000,yt,u_particles(:,:,:,s),u_res(:,:,s),0,0);
%     [eps_ctraj_prop(:,:,s)]=ancestor_trace_LGSS(particles_eps_temp(:,:,:,s),w_temp(:,:,s),indx_temp(:,:,s));
% end
% 
% for s=1:ns
%     eps_prop(:,:,s) = (reshape(eps_ctraj_prop(s,:,:),T,num_prop))';
%     
% end
% 
% for j=1:T
%     eps_collect=[];
%     for k=1:ne
%         eps_collect = [eps_collect,eps_prop(:,j,k)];
%     end
%     covmat_PF(:,:,j) = cov(eps_collect);
%     temp_var = cov(eps_collect);
%     temp_mean = mean(eps_collect);
%     temp_w = 1;
%     
%     temp_w=temp_w.*w_rest;
%     w_mix(:,j)=[w_prior;temp_w];
%     mean_mix(:,:,j)=[zeros(1,ne);temp_mean];
%     covmat_mix(:,:,:,j)=cat(3,eye(ne),temp_var(:,:,1:num_mix));
%     
%     
% end
% 
% %generating initial proposal
% num_prop=100;
% N_prop=1000;
% for s=1:num_prop
%     for t=1:T
%        unif_mix(1,1:N)=sort(rand(1,N));
%        for k=1:num_mix+1
%            id=(unif_mix(1,1:N)>=sum(w_mix(1:(k-1),t))) & (unif_mix(1,1:N)<=sum(w_mix(1:k,t)));
%            id_mix(t,id,s)=k;
%            
%        end
%     end
%     
%    
%     
%     u_particles(:,1:N,:,s)=randn(ne,N,T);
%     u_res(:,1:N,s)=rand(T,N);
%     
%     
%     
% end
% 
% parfor s=1:S
%     [particles_eps(:,:,:,s),w(:,:,s),indx(:,:,s),lik_indv(s,1)]=mix_LGSS_PMMH(B,D,A,C,N,yt,u_particles(:,:,:,s),u_res(:,:,s),...
%         w_mix,mean_mix,covmat_mix,...
%         id_mix(:,:,s),corr_index,disturbance_index);
%     [eps_ctraj(:,:,s)]=ancestor_trace_LGSS(particles_eps(:,:,:,s),w(:,:,s),indx(:,:,s));
% end
% maxw=max(lik_indv);
% var_lik=(var(lik_indv))/S;
% lik_PF=log(mean(exp(lik_indv-maxw)))+maxw; % average of the log likelihood  
% prior_PF=sum(log(unifpdf(vech_A,-1,1))); 
% post_PF=prior_PF+lik_PF;
% 
% D1=num_param;
% V1=0.0001*eye(D1);
% % %V1=covmat;
% scale=0.1;
% accept=0;
% A_rand=0;
% init_time=100;
% % 
% r1_PF=0;
% lik_PF_star=0;
% post_PF_star=0;
% target_accept=0.15;
% for i=1:nloop
% %     
%       %[i,accept,r1_PF]  
% 
%       
%       vech_A = vech(A);
%       R1=mvnrnd(vech_A,scale.*V1);
%       A_star = vech_inv(R1(1,1:end)',dim_states);
%       vech_A_star = vech(A_star);
%       
%       u_particles_star=u_particles;
%       u_res_star=u_res;
%       cindx=randsample(S,1); % select one block to refresh the random numbers
%       
%       
%       u_particles_star(:,:,:,cindx)=corr_param*u_particles(:,:,:,cindx)+sqrt(1-corr_param^2)*randn(dim_states,N,T);
%       u_res_norm=norminv(u_res(:,:,cindx));
%       u_res_norm_star=corr_param*u_res_norm+sqrt(1-corr_param^2)*randn(T,N);
%       u_res_star(:,:,cindx)=normcdf(u_res_norm_star);
%           
%           
% 
%                 
%      parfor s=1:S
%        [particles_eps_star(:,:,:,s),w_star(:,:,s),indx_star(:,:,s),lik_indv_star(s,1)]=mix_LGSS_PMMH(B,D,A_star,C,N,yt,u_particles_star(:,:,:,s),u_res_star(:,:,s),...
%             w_mix,mean_mix,covmat_mix,...
%             id_mix(:,:,s),corr_index,disturbance_index);
%      end
%      if sum(isnan(lik_indv_star))>0 
%         A=A;
%         post_PF=post_PF;
%         u_particles=u_particles;
%         u_res=u_res;
%         particles_eps = particles_eps;
%         w=w;
%         indx=indx;
%      else
%      
%      maxw_star=max(lik_indv_star);
%      var_lik=(var(lik_indv_star))/S;
%      lik_PF_star=log(mean(exp(lik_indv_star-maxw_star)))+maxw_star; % average of the log likelihood
%      
%      prior_PF_star=sum(log(unifpdf(vech_A_star,-1,1))); 
%      
%      post_PF_star=prior_PF_star+lik_PF_star;
%      r1_PF = exp(post_PF_star-post_PF);
%      
%      C1 = min(1,r1_PF);
%      rng('shuffle')
%      A_rand = rand;
%      if A_rand<=C1
%          A=A_star;
%          post_PF=post_PF_star;
%          u_particles=u_particles_star;
%          u_res=u_res_star;
%          accept=accept+1;
%          particles_eps=particles_eps_star;
%          w=w_star;
%          indx=indx_star;
%      end
%      end
% 
%      for s=1:S
%             [eps_ctraj(:,:,s)]=ancestor_trace_LGSS(particles_eps(:,:,:,s),w(:,:,s),indx(:,:,s));
%     
%      end
%      
%      
%      for s=1:ns
%          eps_prop(:,:,s) = (reshape(eps_ctraj(s,:,:),T,num_prop))';
%     
%      end
% 
%     for j=1:T
%         eps_collect=[];
%         for k=1:ne
%             eps_collect = [eps_collect,eps_prop(:,j,k)];
%         end
%         covmat_PF(:,:,j) = cov(eps_collect);
%         temp_var = cov(eps_collect);
%         temp_mean = mean(eps_collect);
%         temp_w = 1;
% 
%         temp_w=temp_w.*w_rest;
%         w_mix(:,j)=[w_prior;temp_w];
%         mean_mix(:,:,j)=[zeros(1,ne);temp_mean];
%         covmat_mix(:,:,:,j)=cat(3,eye(ne),temp_var(:,:,1:num_mix));
%     
%     
%     end
%      
%     thetasave(i,:)=[vech(A)];
% 
%     if i>=init_time
%         V1=cov(thetasave(1:i,:));
%         V1=jitChol(V1);
%         %scale=update_sigma(scale,C1,target_accept,i,D1);
%     end
% %     
%      
% %     for s=1:S
% %        [eps_ctraj(:,:,s)]=ancestor_trace_LGSS(particles_eps(:,:,:,s),w(:,:,s),indx(:,:,s));
% %     end 
%      
%       
%      Post.A(i,:)=vech_A;
%      Post.var(i,1) = var_lik;
%      if mod(i,250)==0
%        %save('/scratch/jz21/dg2271/PMMH_SmallScale_first_order_HS_sysmat_N4.mat','Post'); 
%        save('ParPMMH_LGSS_N100_dim2.mat','Post'); 
%      
%      end
% % 
% % 
% % 
%      
% end

% for i=1:num_states
%     for j=1:i
%         Phi(i,j)=theta^(abs(i-j)+1);
%     end
% end



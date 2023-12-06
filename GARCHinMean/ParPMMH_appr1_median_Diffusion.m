
load('DiffusionDataDaily.mat');
%load('covmat_dim20.mat');
y=y(1:20,:);
T=length(y');
dim_y=size(y,1);
S=101;
N=100;
burn=1000;
nit=29000;
nloop=burn+nit;
%parpool(48)

delta=1;
M=3;
step=delta/M;
indx = obtain_index(dim_y,0);
a = 0.72;
mu=2.36;
tau=0.53;
psi=0.01;


for i=1:dim_y
    for j=1:i
        A(i,j)=psi^(abs(i-j)+1);
    end
end


%vector_A = (A(sub2ind(size(A),[indx(:,1)'],[indx(:,2)'])))';



prior_param.hp_sig2=1;
prior_param.v0=10;
prior_param.s0=1;
prior_param_mean_mu=0;
prior_param_var_mu=10;


ne=dim_y;
ns=ne;
%parpool(28)
corr_param=0.9;
num_param=4;
corr_index=1;
disturbance_index=0;



for s=1:S    
    u_particles(:,1:N,:,s)=randn(ne,N,M*(T-1)+1);
    u_res(:,1:N,s)=rand(M*(T-1)+1,N);
end
parfor s=1:S
    %s
    [~,~,~,lik_indv(:,s)]=smc_diffusion_bootstrap(y,T,M,N,delta,a,tau,mu,A,u_particles(:,:,:,s),u_res(:,:,s),corr_index,disturbance_index);    
end
lik_indv_appr1 = sum(lik_indv);
lik_indv_appr1_sort = sort(lik_indv_appr1);
lik_indv_appr1_median = lik_indv_appr1_sort((S+1)/2);
lik_appr1 = lik_indv_appr1_median;

%end


prior_PF= log(normpdf(psi,0,1))+log(gampdf(a,1,1))+log(gampdf(tau,0.5,0.5))+log(gampdf(mu,1,1));
post_PF=prior_PF+lik_appr1;
jac_PF=log(1/tau)+log(1/a)+log(1/mu);

D1=num_param;
V1=0.001*eye(D1);
%V1=covmat_param;
scale=0.25;
accept=0;
A_rand=0;
target_accept=0.20;
for i=1:nloop
    i
    %accept
    %mu
    %Likelihood at the current value -- lik
    %What we have in the SV code: RW proposal using a transformation
    %vector_A = (A(sub2ind(size(A),[indx(:,1)'],[indx(:,2)'])))';
    theta=[(psi),log(tau),log(a),log(mu)];
    R1=mvnrnd(theta,scale.*V1);
    psi_star = (R1(1,1));
    tau_star = exp(R1(1,2));
    a_star = exp(R1(1,3));
    mu_star = exp(R1(1,4));
    for k=1:dim_y
        for j=1:k
            A_star(k,j)=psi_star^(abs(k-j)+1);
        end
    end
    
    u_particles_star=u_particles;
    u_res_star=u_res;
    cindx=randsample(S,1); % select one block to refresh the random numbers
        
    u_particles_star(:,:,:,cindx)=corr_param*u_particles(:,:,:,cindx)+sqrt(1-corr_param^2)*randn(ne,N,M*(T-1)+1);
    u_res_norm=norminv(u_res(:,:,cindx));
    u_res_norm_star=corr_param*u_res_norm+sqrt(1-corr_param^2)*rand(M*(T-1)+1,N);
    u_res_star(:,:,cindx)=normcdf(u_res_norm_star);
                       
    parfor s=1:S
        [~,~,~,lik_indv_star(:,s)]=smc_diffusion_bootstrap(y,T,M,N,delta,a_star,tau_star,mu_star,A_star,u_particles_star(:,:,:,s),u_res_star(:,:,s),corr_index,disturbance_index);
    end
    
    

    if sum(isnan(sum(lik_indv_star)))>0 
        %sum(sum(sum(sum(isnan(particles_eps_star),4),3),2))>0 
       
    else
    
    lik_indv_appr1_star = sum(lik_indv_star);
    lik_indv_appr1_sort_star = sort(lik_indv_appr1_star);
    lik_indv_appr1_median_star = lik_indv_appr1_sort_star((S+1)/2);
    lik_appr1_star = lik_indv_appr1_median_star;
    
    prior_PF_star= log(normpdf(psi_star))+log(gampdf(a_star,1,1))+log(gampdf(tau_star,0.5,0.5))+log(gampdf(mu_star,1,1));
    jac_PF_star=log(1/tau_star)+log(1/a_star)+log(1/mu_star);
          
    post_PF_star=prior_PF_star+lik_appr1_star;
    r1_PF = exp(post_PF_star-post_PF+jac_PF-jac_PF_star);
    
    C1 = min(1,r1_PF);
    rng('shuffle')
    A_rand = rand;
    if A_rand<=C1
        A = A_star;
        psi=psi_star;
        a = a_star;
        mu = mu_star;
        tau = tau_star;
        accept=accept+1;
        post_PF=post_PF_star;
        jac_PF=jac_PF_star;
        u_particles=u_particles_star;
        u_res=u_res_star;
     end
    
    
   
        
        
        
        
        
    end
    
    
        
    theta=[(psi),log(tau),log(a),log(mu)];
    thetasave(i,:)=theta;
 
    if i>=burn
       %V1=cov(thetasave(1:i,:));
       %V1=jitChol(V1);
       %scale=update_sigma(scale,C1,target_accept,i,6); 
    end
    if scale<0.01
       scale=0.01;
    end
    Post.a(i,1) = a;
    Post.mu(i,1) = mu;
    Post.tau(i,1) = tau;
    Post.psi(i,:) = psi;
        
    %toc
    %delete('PARAMSTORE.mat');
    if mod(i,250)==0
       save(['ParPMMH_appr1_median','dim',num2str(dim_y),'T',num2str(T),'N',num2str(N),'S',num2str(S),'Diffusion_v2.mat'],'Post');
 
    end



end
save(['ParPMMH_appr1_median','dim',num2str(dim_y),'T',num2str(T),'N',num2str(N),'S',num2str(S),'Diffusion_v2.mat'],'Post');
 
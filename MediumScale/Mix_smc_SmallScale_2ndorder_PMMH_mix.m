function [particles_eps,particles_first,particles_second,w,indx,lik]=Mix_smc_SmallScale_2ndorder_PMMH_mix(A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,N,yt,u_particles,u_res,w_mix,mean_mix,covmat_mix,...
    index_mixture,state_var_positions,corr_index,disturbance_index)

ne=size(S2,1);
[~, ns] = size(B);
T         = size(yt,2);
sqrtS2    = R*chol(S2,'lower');
num_mix = size(w_mix,1);

particles_first=zeros(ns,N,T);
particles_second=zeros(ns,N,T);
particles_eps=zeros(ne,N,T);
w=zeros(T,N);
indx=zeros(T,N);
%w1_mix=0.15;
%w2_mix=0.80;
%w3_mix=1-w1_mix-w2_mix;

x_f_init=zeros(size(Phi,2),1);  %Initilization
x_s_init=zeros(size(Phi,2),1);
t=1;
yy=yt(:,t);
particles_eps_temp=[];
for k=1:num_mix
    n_mix(k,1)=sum(index_mixture(t,1:N)==k);
end
for k=1:num_mix
    if k==1
        temp=mean_mix(k,:,t)'+chol(covmat_mix(:,:,k,t),'lower')*u_particles(:,1:n_mix(k,1),t);
    else
        temp=mean_mix(k,:,t)'+chol(covmat_mix(:,:,k,t),'lower')*u_particles(:,sum(n_mix(1:k-1,1))+1:sum(n_mix(1:k,1)),t);
    end
    particles_eps_temp=[particles_eps_temp,temp];
end

particles_eps(:,:,t)=particles_eps_temp;
% n1=sum(index_mixture(t,1:N)==1);
% n2=sum(index_mixture(t,1:N)==2);
% n3=sum(index_mixture(t,1:N)==3);
% particles_eps_temp1=eps_ctraj(:,t)+chol(covmat_eps(:,:,t),'lower')*u_particles(:,1:n1,t);
% particles_eps_temp2=mean_eps(t,:)'+chol(covmat_eps(:,:,t),'lower')*u_particles(:,n1+1:n1+n2,t);
% particles_eps_temp3=chol(eye(ne),'lower')*u_particles(:,n1+n2+1:n1+n2+n3,t);
% particles_eps(:,:,t)=[particles_eps_temp1,particles_eps_temp2,particles_eps_temp3];
particles_first(:,:,t)=Phi*x_f_init + sqrtS2*particles_eps(:,:,t);
x_f(:,:,t)=particles_first(state_var_positions,:,t);
particles_second(:,:,t)=0.5*omega_ss_delta+Phi*x_s_init+...
    sqrtS2*particles_eps(:,1:N,t)+0.5*Phi11*(kron(x_f_init,ones(size(Phi,2),1)).*repmat(x_f_init,size(Phi,2),1))+...
    Phi12*(kron(x_f_init,ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),size(Phi,2),1))+...
    0.5*R00*(kron(chol(S2,'lower')*particles_eps(:,1:N,t),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),ne,1));
x_s(:,:,t)=particles_second(state_var_positions,:,t);
perror = repmat(yy-A,1,N)-B*particles_second(:,:,t);
logw1 = logmvnpdf(perror',zeros(1,size(yy,2)),H);
logw2 = logmvnpdf(particles_eps(:,:,t)',zeros(1,size(yy,2)),eye(ne));

for k=1:num_mix
    temp_dens(:,k)=w_mix(k,t).*mvnpdf(particles_eps(:,:,t)',mean_mix(k,:,t),covmat_mix(:,:,k,t));
end
logw3=log(sum(temp_dens,2));
%logw3 = log(w1_mix.*mvnpdf(particles_eps(:,:,t)',eps_ctraj(:,t)',covmat_eps(:,:,t))+...
%    w2_mix.*mvnpdf(particles_eps(:,:,t)',mean_eps(t,:),covmat_eps(:,:,t))+...
%    w3_mix.*mvnpdf(particles_eps(:,:,t)',zeros(1,ne),eye(ne)));
logw=logw1+logw2-logw3';
w(t,:)=exp(logw-max(logw));
lik(t,1)=log(mean(w(t,:)))+max(logw);
w(t,:)=w(t,:)./sum(w(t,:));  
 
for t=2:T
    yy=yt(:,t);
    %particles_com=[x_f(:,1:N,t-1);x_s(:,1:N,t-1)];
    %particles_com=[particles_first(:,1:N,t-1);particles_second(:,1:N,t-1)];
    %particles_com=[particles_second(:,1:N,t-1)];
    
    if corr_index==1
       if disturbance_index==1    
          particles_com=[particles_eps(:,1:N,t-1)];
       else
          particles_com=[particles_first(:,1:N,t-1);particles_second(:,1:N,t-1)];
       end 
       %[indx_temp(t,:)]=rs_multinomial_corr2_PMMH_2ndorder_anti(particles_com,w(t-1,:),u_res(t,:));
       [indx_temp(t,:)]=rs_multinomial_corr2_PMMH_2ndorder_simplified(particles_com,w(t-1,:),u_res(t,:));
    
    else
       [indx_temp(t,:)]=rs_multinomial(w(t-1,:));
    end
    indx(t,:)=[indx_temp(t,:)];
    
    %indx(t,:)=[indx_temp(t,:),indx_temp(t,:)];
    particles_eps_temp=[];
    for k=1:num_mix
        n_mix(k,1)=sum(index_mixture(t,1:N)==k);
    end
    for k=1:num_mix
        if k==1
            temp=mean_mix(k,:,t)'+chol(covmat_mix(:,:,k,t),'lower')*u_particles(:,1:n_mix(k,1),t);
        else
            temp=mean_mix(k,:,t)'+chol(covmat_mix(:,:,k,t),'lower')*u_particles(:,sum(n_mix(1:k-1,1))+1:sum(n_mix(1:k,1)),t);
        end
        particles_eps_temp=[particles_eps_temp,temp];
    end
    particles_eps(:,:,t)=particles_eps_temp;

%     n1=sum(index_mixture(t,1:N)==1);
%     n2=sum(index_mixture(t,1:N)==2);
%     n3=sum(index_mixture(t,1:N)==3); 
%      
%     particles_eps_temp1=eps_ctraj(:,t)+chol(covmat_eps(:,:,t),'lower')*u_particles(:,1:n1,t);
%     particles_eps_temp2=mean_eps(t,:)'+chol(covmat_eps(:,:,t),'lower')*u_particles(:,n1+1:n1+n2,t);
%     particles_eps_temp3=chol(eye(ne),'lower')*u_particles(:,n1+n2+1:n1+n2+n3,t);
%     particles_eps(:,:,t)=[particles_eps_temp1,particles_eps_temp2,particles_eps_temp3];
    particles_first(:,:,t)=Phi*x_f(:,indx(t,1:N),t-1)+sqrtS2*particles_eps(:,1:N,t);
    x_f(:,:,t)=particles_first(state_var_positions,:,t);
    particles_second(:,1:N,t)=0.5*omega_ss_delta+Phi*x_s(:,indx(t,1:N),t-1)+...
        sqrtS2*particles_eps(:,1:N,t)+0.5*Phi11*(kron(x_f(:,indx(t,1:N),t-1),ones(size(Phi,2),1)).*repmat(x_f(:,indx(t,1:N),t-1),size(Phi,2),1))+...
        Phi12*(kron(x_f(:,indx(t,1:N),t-1),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),size(Phi,2),1))+...
        0.5*R00*(kron(chol(S2,'lower')*particles_eps(:,1:N,t),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),ne,1));
    x_s(:,:,t)=particles_second(state_var_positions,:,t);
    perror = repmat(yy-A,1,N)-B*particles_second(:,:,t);
    logw1=logmvnpdf(perror',zeros(1,size(yy,2)),H);
    logw2=logmvnpdf(particles_eps(:,:,t)',zeros(1,size(yy,2)),eye(ne));
    
    for k=1:num_mix
        temp_dens(:,k)=w_mix(k,t).*mvnpdf(particles_eps(:,:,t)',mean_mix(k,:,t),covmat_mix(:,:,k,t));
    end
    logw3=log(sum(temp_dens,2));
    
%     logw3 = log(w1_mix.*mvnpdf(particles_eps(:,:,t)',eps_ctraj(:,t)',covmat_eps(:,:,t))+...
%           w2_mix.*mvnpdf(particles_eps(:,:,t)',mean_eps(t,:),covmat_eps(:,:,t))+...
%           w3_mix.*mvnpdf(particles_eps(:,:,t)',zeros(1,ne),eye(ne)));
    logw = logw1+logw2-logw3';
    w(t,:)=exp(logw-max(logw));
    lik(t,1)=log(mean(w(t,:)))+max(logw);
    w(t,:)=w(t,:)./sum(w(t,:)); 

end
end



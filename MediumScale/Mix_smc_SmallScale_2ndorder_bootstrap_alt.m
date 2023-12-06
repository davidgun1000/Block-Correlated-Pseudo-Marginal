function [particles_eps,particles_first,particles_second,w,indx,lik]=Mix_smc_SmallScale_2ndorder_bootstrap_alt(A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,N,yt,state_var_positions)

ne=size(S2,1);
[~, ns] = size(B);
T         = size(yt,2);
sqrtS2    = R*chol(S2,'lower');

particles_first=zeros(ns,N,T);
particles_second=zeros(ns,N,T);
particles_eps=zeros(ne,N,T);
w=zeros(T,N);
indx=zeros(T,N);
x_f_init=zeros(size(Phi,2),1);  %Initilization
x_s_init=zeros(size(Phi,2),1);
t=1;
yy=yt(:,t);
particles_eps(:,:,t)=randn(ne,N);
particles_first(:,:,t)=Phi*x_f_init + sqrtS2*particles_eps(:,:,t);
x_f(:,:,t)=particles_first(state_var_positions,:,t);
particles_second(:,:,t)=0.5*omega_ss_delta+Phi*x_s_init+...
    sqrtS2*particles_eps(:,1:N,t)+0.5*Phi11*(kron(x_f_init,ones(size(Phi,2),1)).*repmat(x_f_init,size(Phi,2),1))+...
    Phi12*(kron(x_f_init,ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),size(Phi,2),1))+...
    0.5*R00*(kron(chol(S2,'lower')*particles_eps(:,1:N,t),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),ne,1));
x_s(:,:,t)=particles_second(state_var_positions,:,t);
perror = repmat(yy-A,1,N)-B*particles_second(:,:,t);
logw1 = logmvnpdf(perror',zeros(1,size(yy,1)),H);
logw = logw1;
w(t,:)=exp(logw-max(logw));
lik(t,1)=log(mean(w(t,:)))+max(logw);
w(t,:)=w(t,:)./sum(w(t,:)); 
 
for t=2:T
    yy=yt(:,t);
    [indx(t,:)]=rs_multinomial(w(t-1,:));
    particles_eps(:,:,t)=randn(ne,N);
    particles_first(:,:,t)=Phi*x_f(:,indx(t,1:N),t-1)+sqrtS2*particles_eps(:,1:N,t);
    x_f(:,:,t)=particles_first(state_var_positions,:,t);
    particles_second(:,1:N,t)=0.5*omega_ss_delta+Phi*x_s(:,indx(t,1:N),t-1)+...
        sqrtS2*particles_eps(:,1:N,t)+0.5*Phi11*(kron(x_f(:,indx(t,1:N),t-1),ones(size(Phi,2),1)).*repmat(x_f(:,indx(t,1:N),t-1),size(Phi,2),1))+...
        Phi12*(kron(x_f(:,indx(t,1:N),t-1),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),size(Phi,2),1))+...
        0.5*R00*(kron(chol(S2,'lower')*particles_eps(:,1:N,t),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),ne,1));
    x_s(:,:,t)=particles_second(state_var_positions,:,t);
    perror = repmat(yy-A,1,N)-B*particles_second(:,:,t);
    logw1=logmvnpdf(perror',zeros(1,size(yy,1)),H);
    logw = logw1;
    w(t,:)=exp(logw-max(logw));
    lik(t,1)=log(mean(w(t,:)))+max(logw);
    w(t,:)=w(t,:)./sum(w(t,:)); 
end
end


% 
% %initialisation
% t=1;
% x0 = zeros(ns,1);
% P0=init_P0;
% P0=jitChol(P0);
% P0=topdm(P0);
% temp_s=x0;
% temp_P=P0;
% particles(:,1:N,t)=repmat(temp_s,1,N)+chol(temp_P)'*u_init;
% 
% t=2;
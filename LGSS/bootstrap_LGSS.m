function [particles,w,indx,lik]=bootstrap_LGSS(y,D,Z,sigma_u,C,G,M,sigma_e,u_particles_init,u_particles,u_res,N)

T=size(y,2);
[~, ns] = size(Z);
particles=zeros(ns,N,T);
w=zeros(T,N);
indx=zeros(T,N);
t=1;
xinit=zeros(ns,1); % initial target state
Vinit=reshape((eye(ns^2)-kron(G,G))\reshape(M*sigma_e*M',ns^2,1),ns,ns); % filter's initial variance
particles_init=xinit+chol(Vinit,'lower')*u_particles_init;

yy=y(:,t);
particles(:,:,t) = G*particles_init + chol(sigma_e,'lower')*u_particles(:,:,t);
perror = repmat(yy,1,N)-particles(:,:,t);
logw=logmvnpdf(perror',zeros(1,size(yy,2)),sigma_u);
w(t,:)=exp(logw-max(logw));
lik(t,1)=log(mean(w(t,:)))+max(logw);
w(t,:)=w(t,:)./sum(w(t,:)); 
 
for t=2:T
    yy=y(:,t);
    particles_com=particles(:,:,t-1);
    [indx(t,:)]=rs_multinomial(w(t-1,:));
    particles(:,:,t) = G*particles(:,indx(t,1:N),t-1) + chol(sigma_e,'lower')*u_particles(:,:,t);
    perror = repmat(yy,1,N)-particles(:,:,t);
    logw=logmvnpdf(perror',zeros(1,size(yy,2)),sigma_u);
    w(t,:)=exp(logw-max(logw));
    lik(t,1)=log(mean(w(t,:)))+max(logw);
    w(t,:)=w(t,:)./sum(w(t,:)); 

     
     
     %     perror=repmat(yy,1,N)-particles(:,:,t);
%     
%     logw=logmvnpdf(perror',zeros(1,size(yy,2)),H);
%     w(t,:)=exp(logw-max(logw));
%     lik=lik+log(mean(w(t,:)))+max(logw);
%     w(t,:)=w(t,:)./sum(w(t,:)); 
end

end
% ne=size(S2,1);
% [~, ns] = size(B);
% T         = size(yt,2);
% sqrtS2    = R*chol(S2,'lower');
% 
% particles_first=zeros(ns,N,T);
% particles_second=zeros(ns,N,T);
% particles_eps=zeros(ne,N,T);
% w=zeros(T,N);
% indx=zeros(T,N);
% x_f_init=zeros(size(Phi,2),1);  %Initilization
% x_s_init=zeros(size(Phi,2),1);
% t=1;
% yy=yt(:,t);
% particles_eps(:,:,t)=randn(ne,N);
% particles_first(:,:,t)=Phi*x_f_init + sqrtS2*particles_eps(:,:,t);
% x_f(:,:,t)=particles_first(state_var_positions,:,t);
% particles_second(:,:,t)=0.5*omega_ss_delta+Phi*x_s_init+...
%     sqrtS2*particles_eps(:,1:N,t)+0.5*Phi11*(kron(x_f_init,ones(size(Phi,2),1)).*repmat(x_f_init,size(Phi,2),1))+...
%     Phi12*(kron(x_f_init,ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),size(Phi,2),1))+...
%     0.5*R00*(kron(chol(S2,'lower')*particles_eps(:,1:N,t),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),ne,1));
% x_s(:,:,t)=particles_second(state_var_positions,:,t);
% perror = repmat(yy,1,N)-B*particles_second(:,:,t);
% logw1 = logmvnpdf(perror',zeros(1,size(yy,2)),H);
% logw = logw1;
% w(t,:)=exp(logw-max(logw));
% lik=log(mean(w(t,:)))+max(logw);
% w(t,:)=w(t,:)./sum(w(t,:)); 
%  
% for t=2:T
%     yy=yt(:,t);
%     [indx(t,:)]=rs_multinomial(w(t-1,:));
%     particles_eps(:,:,t)=randn(ne,N);
%     particles_first(:,:,t)=Phi*x_f(:,indx(t,1:N),t-1)+sqrtS2*particles_eps(:,1:N,t);
%     x_f(:,:,t)=particles_first(state_var_positions,:,t);
%     particles_second(:,1:N,t)=0.5*omega_ss_delta+Phi*x_s(:,indx(t,1:N),t-1)+...
%         sqrtS2*particles_eps(:,1:N,t)+0.5*Phi11*(kron(x_f(:,indx(t,1:N),t-1),ones(size(Phi,2),1)).*repmat(x_f(:,indx(t,1:N),t-1),size(Phi,2),1))+...
%         Phi12*(kron(x_f(:,indx(t,1:N),t-1),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),size(Phi,2),1))+...
%         0.5*R00*(kron(chol(S2,'lower')*particles_eps(:,1:N,t),ones(ne,1)).*repmat(chol(S2,'lower')*particles_eps(:,1:N,t),ne,1));
%     x_s(:,:,t)=particles_second(state_var_positions,:,t);
%     perror = repmat(yy,1,N)-B*particles_second(:,:,t);
%     logw1=logmvnpdf(perror',zeros(1,size(yy,2)),H);
%     logw = logw1;
%     w(t,:)=exp(logw-max(logw));
%     lik=lik+log(mean(w(t,:)))+max(logw);
%     w(t,:)=w(t,:)./sum(w(t,:)); 
% end



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
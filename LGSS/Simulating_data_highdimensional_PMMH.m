%%Simulating data

T=300;
dim_states=15;
sigma_u = eye(dim_states);
sigma_e = eye(dim_states);
Z = eye(dim_states);
M=eye(dim_states);
%temp_A = rand(dim_states*(dim_states+1)/2,1);
%A = vech_inv(temp_A,dim_states);
%A(1,:) = [0.9,0];
%A(2,:) = [0.3,0.7];
theta=0.4;
for i=1:dim_states
    for j=1:i
        G(i,j)=theta^(abs(i-j)+1);
    end
end
x_init_mean=zeros(dim_states,1); % initial target state
V_init=reshape((eye(dim_states^2)-kron(G,G))\reshape(M*sigma_e*M',dim_states^2,1),dim_states,dim_states); % filter's initial variance



t=1;
x_init = x_init_mean+chol(V_init,'lower')*randn(dim_states,1);
xt(:,t) = G*x_init + chol(sigma_e,'lower')*randn(dim_states,1);
yt(:,t) = Z*xt(:,t) + chol(sigma_u,'lower')*randn(dim_states,1);

for t=2:T
    eps=randn(dim_states,1);
    xt(:,t)=G*xt(:,t-1)+chol(sigma_e,'lower')*randn(dim_states,1);
    yt(:,t)=Z*xt(:,t)+chol(sigma_u,'lower')*randn(dim_states,1);
end
xt_true = xt;
G_true = G;
save(['simdata',num2str(dim_states),'_',num2str(T),'.mat'],'yt','xt_true','G_true','Z','sigma_u','sigma_e','dim_states','T');

%     %yt(:,t)=B*xt(:,t)+chol(H,'lower')*randn(num_observed,1);    
%     yt(:,t)=xt(:,t)+chol(H,'lower')*randn(num_states,1);    
% 


% T=100;
% tau_true=1;
% sig_true=1;
% num_states=10;
% num_observed=2;
% theta=0.2;
% 
% for i=1:num_states
%     for j=1:i
%         phi_true(i,j)=theta^(abs(i-j)+1);
%     end
% end
% t=1;
% B=zeros(num_observed,num_states);
% B(1,1)=1;
% B(2,2)=1;
% H=0.5*eye(num_states);
% %H=diag([sig_true^2*ones(num_states,1)]);
% %xt(:,t)=sqrt(tau_true/(1-phi_true^2)).*randn(num_states,1);
% xt(:,t)=randn(num_states,1);
% 
% %yt(:,t)=B*xt(:,t)+chol(H,'lower')*randn(num_observed,1);
% yt(:,t)=xt(:,t)+chol(H,'lower')*randn(num_states,1);
% 
% for t=2:T
%     eps=randn(num_states,1);
%     xt(:,t)=phi_true*xt(:,t-1)+tau_true.*eps;
%     %yt(:,t)=B*xt(:,t)+chol(H,'lower')*randn(num_observed,1);    
%     yt(:,t)=xt(:,t)+chol(H,'lower')*randn(num_states,1);    
% 
% end
% 


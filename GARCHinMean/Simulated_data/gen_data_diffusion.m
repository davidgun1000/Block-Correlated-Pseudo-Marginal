%generate data
data=6;
delta=1;
M=3;
step=delta/M;
dim_y=20;
T=100;
alpha_true=2;
mu_true=1.81;
tau_true=0.38;
psi_true=0.01;
t=1;

%[indx]=obtain_index(dim_y,0);

h_true(:,t)=randn(dim_y,1);
%temp = unifrnd(-0.2,0.2,1,1);
%A_true=zeros(dim_y);
%A_true(sub2ind(size(A_true),indx(:,1),indx(:,2))) = temp';

for i=1:dim_y
    for j=1:i
        A_true(i,j)=psi_true^(abs(i-j)+1);
    end
end


%A_true = vech_inv(vech_A,dim_y);
sigma_true = diag([exp(h_true(:,t))]);
r=1;
y(:,r) = A_true*h_true(:,t) + chol(sigma_true,'lower')*randn(dim_y,1);

for t=2:(M*(T-1)+1)
    for i=1:dim_y
       h_true(i,t) = h_true(i,t-1) + step.*(alpha_true.*(mu_true-exp(h_true(i,t-1))).*exp(-h_true(i,t-1))-(tau_true/2))+sqrt(tau_true)*sqrt(step)*randn;         
    end

%     sigma_true = diag([exp(h_true(:,t))]);
%     y(:,t) = A_true*h_true(:,t)  + chol(sigma_true,'lower')*randn(dim_y,1);



    %h_true(1,t) = h_true(1,t-1) + step.*(alpha_true.*(mu_true-exp(h_true(1,t-1))).*exp(-h_true(1,t-1))-(tau_true/2))+sqrt(tau_true)*sqrt(step)*randn;
    %h_true(2,t) = h_true(2,t-1) + step.*(alpha_true.*(mu_true-exp(h_true(2,t-1))).*exp(-h_true(2,t-1))-(tau_true/2))+sqrt(tau_true)*sqrt(step)*randn;
    
    if mod(t,M)==1
       r=r+1;
       sigma_true = diag([exp(h_true(:,t))]);
       y(:,r) = A_true*h_true(:,t)  + chol(sigma_true,'lower')*randn(dim_y,1);
        
        
    end
    
end

save(['sim_GARCH_diffusion',num2str(dim_y),'_','data',num2str(data),'new','.mat'],'alpha_true','mu_true','tau_true','A_true','h_true','y','psi_true');

%vech_A = (A(sub2ind(size(A),[indx(:,1)'],[indx(:,2)'])))';
%temp_mat(sub2ind(size(temp_mat),indx(:,1),indx(:,2)))=temp_vec';
%temp_vec=grad_T(sub2ind(size(grad_T),[indx(:,1)'],[indx(:,2)']));  
%vech_A=unifrnd(-0.2,0.2,length(indx),1);

%A0_11_true=0.1;
%A0_22_true=0.1;
%A0_21_true=0.05;
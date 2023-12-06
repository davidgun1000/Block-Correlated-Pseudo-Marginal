function [s_fore_first_new, s_fore_second_new, eps_new, ind_acpt] = mutation_RWMH_multistep(s0_fore_first, s0_fore_second, ...
    s0_up_first,s0_up_second, eps0, A,B,H, Phi, R, S2, Phi11, Phi12, R00, omega_ss_delta, yy, tune)
% RWMH for mutation step
%mutation_RWMH_multistep(s_fore_first(:,nn), s_fore_second(:,nn),...
%                     s_up_first(:,nn),s_up_second(:,nn), eps_up(:,nn), B, H/phi_new, Phi, R, S2, Phi11, Phi12, R00, omega_ss_delta, yy, tune);
ind_acpt = 0;
ne=size(S2,1);
[~,ns]=size(B);
sqrtS2=R*chol(S2,'lower');
state_var_positions=[2;3;4;5];
% iterate over MH steps
for nn = 1:tune.N_mh

    % RW proposal
    epsx = eps0 + tune.c*chol(tune.R)'*randn(size(tune.R,1),1);
    sx_fore_first   = Phi*s0_up_first + sqrtS2*epsx;
    sx_fore_second  = 0.5*omega_ss_delta+Phi*s0_up_second+... 
       sqrtS2*epsx+0.5*Phi11*(kron(s0_up_first,ones(size(Phi,2),1)).*repmat(s0_up_first,size(Phi,2),1))+... 
       Phi12*(kron(s0_up_first,ones(ne,1)).*repmat(chol(S2,'lower')*epsx,size(Phi,2),1))+...
        0.5*R00*(kron(chol(S2,'lower')*epsx,ones(ne,1)).*repmat(chol(S2,'lower')*epsx,ne,1)); 
    
     perrorx = yy -A - B*sx_fore_second;
     perror0 = yy -A - B*s0_fore_second;
    
     postx = log(mvnpdf(perrorx', zeros(1,size(yy,2)), H)*mvnpdf(epsx'));
     post0 = log(mvnpdf(perror0', zeros(1,size(yy,2)), H)*mvnpdf(eps0'));

    % Accept/Reject
     alp = exp(postx - post0); % this is RW, so q is canceled out
     if rand < alp % accept
        eps_new    = epsx;
        s_fore_first_new = sx_fore_first;
        s_fore_second_new= sx_fore_second;
        ind_acpt   = ind_acpt + 1;
     else % reject
        eps_new    = eps0;
        s_fore_first_new = s0_fore_first;
        s_fore_second_new= s0_fore_second;
        ind_acpt   = 0;
     end
    
     eps0 = eps_new;
     s0_fore_first=s_fore_first_new;
     s0_fore_second=s_fore_second_new;
     
     %eps0 = eps_new;
     %s0_fore = s_fore_new;
     
     %s0_up_first = s_fore_first_new(state_var_positions,:);
     %s0_up_second= s_fore_second_new(state_var_positions,:);   
end

ind_acpt = ind_acpt/tune.N_mh;
end

% % outside of function
% parasim(i,j,:) = ind_para;
% loglh(j)       = ind_loglh;
% temp_acpt(j,1) = ind_acpt;

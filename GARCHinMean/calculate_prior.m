function [prior_param]=calculate_prior(param_col,length_A,length_mu,length_tau,length_a,prior_param)
   
   vector_A = param_col(1,1:length_A)';
   mu = exp(param_col(1,length_A+1:length_A+length_mu))';
   tau = exp(param_col(1,length_A+length_mu+1:length_A+length_mu+length_tau)');
   a = exp(param_col(1,length_A+length_mu+length_tau+1:length_A+length_mu+length_tau+length_a)');

   prior_param = sum(log(normpdf(vector_A)))+sum(log_IG_PDF_used(a,prior_param.v0/2,prior_param.s0/2))+...
       sum(log_IG_PDF_used(tau,prior_param.v0/2,prior_param.s0/2))+sum(log(normpdf(mu,prior_param.mean_mu,prior_param.std_mu)));
    

end
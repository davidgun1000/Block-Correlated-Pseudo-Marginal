function [jacobian_param]=calculate_jacobian(param_col,length_A,length_mu,length_tau,length_a)
   
   vector_A = param_col(1,1:length_A)';
   mu = exp(param_col(1,length_A+1:length_A+length_mu))';
   tau = exp(param_col(1,length_A+length_mu+1:length_A+length_mu+length_tau)');
   a = exp(param_col(1,length_A+length_mu+length_tau+1:length_A+length_mu+length_tau+length_a)');

   jacobian_param = sum(log(1./tau))+sum(log(1./a))+sum(log(1./mu));
   
       

end
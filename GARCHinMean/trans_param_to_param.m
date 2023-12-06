function [param] = trans_param_to_param(trans_param,length_A,length_mu,length_tau,length_a)

    length_param = length(trans_param);
    
    trans_param = trans_param';
    ind_select = [1:1:length_param]';

    for i=1:length_param
        
        if ind_select(i,1)<=length_A
           param(i,1) = trans_param(i,1); 
        end
        
        if ind_select(i,1)>length_A & ind_select(i,1)<=length_A+length_mu
           param(i,1) = exp(trans_param(i,1));
        end

        if ind_select(i,1)>length_A+length_mu & ind_select(i,1)<=length_A+length_mu+length_tau
           param(i,1) = exp(trans_param(i,1));
        end

        if ind_select(i,1)>length_A+length_mu+length_tau & ind_select(i,1)<=length_A+length_mu+length_tau+length_a
           param(i,1) = exp(trans_param(i,1));
        end
        
    end

end

% switch ind_select(i,1)
%         
%             case 1
%                temp = exp(trans_param(i,1));            
%             case 2
%                temp = exp(trans_param(i,1));
%             case 3
%                temp = exp(trans_param(i,1)); 
%             case 4 
%                temp = logitcdf(trans_param(i,1));
%             case 5
%                temp = logitcdf(trans_param(i,1));
%             case 6
%                temp = logitcdf(trans_param(i,1));
%             case 7
%                temp = exp(trans_param(i,1)); 
%             case 8 
%                temp = exp(trans_param(i,1)); 
%             case 9
%                temp = exp(trans_param(i,1)); 
%             case 10
%                temp = exp(trans_param(i,1));  
%             case 11
%                temp = exp(trans_param(i,1));
%             case 12
%                temp = trans_param(i,1);
%             case 13
%                 temp = exp(trans_param(i,1)); 
%             case 14
%                 temp = exp(trans_param(i,1)); 
%             case 15
%                  temp = exp(trans_param(i,1)); 
%         end
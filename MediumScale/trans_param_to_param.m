function [param] = trans_param_to_param(trans_param,ind_select)

    length_param = length(trans_param);
    trans_param = trans_param';
    for i=1:length_param
        
        switch ind_select(i,1)
        
            case 1
               temp = logitcdf(trans_param(i,1));            
            case 2
               temp = logitcdf(trans_param(i,1));
            case 3
               temp = logitcdf(trans_param(i,1)); 
            case 4 
               temp = exp(trans_param(i,1));
            case 5
               temp = exp(trans_param(i,1));
            case 6
               temp = exp(trans_param(i,1));
            case 7
               temp = exp(trans_param(i,1)); 
            case 8 
               temp = exp(trans_param(i,1)); 
            case 9
               temp = exp(trans_param(i,1)); 
            case 10
               temp = exp(trans_param(i,1));  
            case 11
               temp = exp(trans_param(i,1));
            case 12
               temp = exp(trans_param(i,1));
            case 13
                temp = exp(trans_param(i,1)); 
            case 14
                temp = (trans_param(i,1)); 
            case 15
                 temp = (trans_param(i,1)); 
            case 16
                temp = (trans_param(i,1));
            case 17
                temp = logitcdf(trans_param(i,1));
            case 18
                temp = exp(trans_param(i,1));
            case 19
                temp = trans_param(i,1);
            case 20
                temp = trans_param(i,1);
            case 21
                temp = exp(trans_param(i,1));
            case 22
                temp = exp(trans_param(i,1));
            case 23
                temp = logitcdf(trans_param(i,1));
            case 24
                temp = logitcdf(trans_param(i,1));
            case 25
                temp = exp(trans_param(i,1));
            case 26
                temp = trans_param(i,1);
            case 27
                temp = trans_param(i,1);
            case 28
                temp = trans_param(i,1);
        end
        
        param(i,1) = temp;
        
    end





end
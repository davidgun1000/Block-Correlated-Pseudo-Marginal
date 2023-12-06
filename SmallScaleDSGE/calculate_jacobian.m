function [log_jac] = calculate_jacobian(param_col)

    log_jac = log(1/param_col(1,1))+log(1/param_col(2,1))+log(1/param_col(3,1))+...
    log((1/param_col(4,1))+(1/(1-param_col(4,1))))+log((1/param_col(5,1))+(1/(1-param_col(5,1))))+log((1/param_col(6,1))+(1/(1-param_col(6,1))))+...
    log(1/param_col(7,1))+log(1/param_col(8,1))+log(1/param_col(9,1))+log(1/param_col(11,1))+log(1/param_col(10,1))+...
    log(1/param_col(13,1))+log(1/param_col(14,1))+log(1/param_col(15,1));



end
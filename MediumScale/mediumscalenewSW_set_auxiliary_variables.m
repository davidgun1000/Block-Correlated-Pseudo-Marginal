function y = mediumscalenewSW_set_auxiliary_variables(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

y(50)=exp(x(4))*(y(12)*exp(x(4))/y(12)-1)*exp(y(20))*y(13)*params(2)*params(9);
y(51)=exp(x(4))*(exp(x(4))*y(38)/y(38)-1)*exp(y(20))*y(37)*params(2)*params(9);

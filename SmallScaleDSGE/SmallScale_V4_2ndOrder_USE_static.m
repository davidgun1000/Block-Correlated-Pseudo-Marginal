function [residual, g1, g2, g3] = SmallScale_V4_2ndOrder_USE_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 8, 1);

%
% Model equations
%

T35 = exp(y(6))*(1-0.5/params(2))+exp(params(4))*0.5/params(2);
T53 = exp(y(5))^(1-params(9));
T54 = exp(y(4))^params(9);
T63 = exp(y(6))/exp(params(4));
T66 = exp(params(4))*params(13)/params(1)*T63^params(6);
T72 = (1-params(2))^(1/params(5))*exp(y(7));
T73 = exp(y(1))/T72;
T75 = T73^params(7);
T82 = params(3)*0.5*exp(params(4))^2;
lhs =params(1)*exp(y(4))/(exp(y(6))*exp(y(8)));
rhs =1;
residual(1)= lhs-rhs;
lhs =params(3)*(exp(y(6))-exp(params(4)))*T35-exp(y(6))*(exp(y(6))-exp(params(4)))*params(1)*params(3)+1/params(2)*(1-exp(y(3))^params(5));
rhs =1;
residual(2)= lhs-rhs;
lhs =exp(y(4));
rhs =T53*T54*exp(x(1));
residual(3)= lhs-rhs;
lhs =exp(y(5));
rhs =T66*T75;
residual(4)= lhs-rhs;
lhs =exp(y(3))/exp(y(1));
rhs =1/exp(y(7))-T82*(T63-1)^2;
residual(5)= lhs-rhs;
lhs =y(7);
rhs =(1-params(8))*params(10)+y(7)*params(8)+x(2);
residual(6)= lhs-rhs;
lhs =y(8);
rhs =(1-params(17))*params(14)+y(8)*params(17)+x(3);
residual(7)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(1));
residual(8)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(8, 8);

  %
  % Jacobian matrix
  %

T109 = getPowerDeriv(T73,params(7),1);
T135 = (-(params(1)*exp(y(4))*exp(y(6))*exp(y(8))))/(exp(y(6))*exp(y(8))*exp(y(6))*exp(y(8)));
  g1(1,4)=params(1)*exp(y(4))/(exp(y(6))*exp(y(8)));
  g1(1,6)=T135;
  g1(1,8)=T135;
  g1(2,3)=1/params(2)*(-(exp(y(3))*getPowerDeriv(exp(y(3)),params(5),1)));
  g1(2,6)=T35*exp(y(6))*params(3)+params(3)*(exp(y(6))-exp(params(4)))*exp(y(6))*(1-0.5/params(2))-(exp(y(6))*(exp(y(6))-exp(params(4)))*params(1)*params(3)+exp(y(6))*exp(y(6))*params(1)*params(3));
  g1(3,4)=exp(y(4))-exp(x(1))*T53*exp(y(4))*getPowerDeriv(exp(y(4)),params(9),1);
  g1(3,5)=(-(exp(x(1))*T54*exp(y(5))*getPowerDeriv(exp(y(5)),1-params(9),1)));
  g1(4,1)=(-(T66*T73*T109));
  g1(4,5)=exp(y(5));
  g1(4,6)=(-(T75*exp(params(4))*params(13)/params(1)*T63*getPowerDeriv(T63,params(6),1)));
  g1(4,7)=(-(T66*T109*(-(exp(y(1))*T72))/(T72*T72)));
  g1(5,1)=(-(exp(y(3))*exp(y(1))))/(exp(y(1))*exp(y(1)));
  g1(5,3)=exp(y(3))/exp(y(1));
  g1(5,6)=T82*T63*2*(T63-1);
  g1(5,7)=(-((-exp(y(7)))/(exp(y(7))*exp(y(7)))));
  g1(6,7)=1-params(8);
  g1(7,8)=1-params(17);
  g1(8,1)=(-exp(y(1)));
  g1(8,2)=exp(y(2));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,64);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,512);
end
end
end
end

function [residual, g1, g2, g3] = medscalenew_static(y, x, params)
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

residual = zeros( 30, 1);

%
% Model equations
%

T10 = 1/y(22);
T16 = (y(2)-y(2)*params(7)/y(22))^(-1);
T38 = y(5)/y(6)-1;
T66 = y(1)*y(10)*params(6)/params(11);
T70 = params(29)*y(10)^params(8);
T75 = T70/y(1)-y(9)*(params(6)-1)/params(6);
T89 = y(12)*exp(x(4))/y(12)-1;
T95 = y(13)*exp(y(20))*(1-y(12)*exp(x(4))*params(9)*T89/y(12));
T103 = T89^2;
T133 = y(10)*y(9)*y(22)*params(3)/(1-params(3));
T144 = params(24)^params(16);
T152 = params(46)*params(24)^params(18);
T167 = y(14)^(-1);
T182 = (y(16)*y(18)/y(22))^params(3);
T183 = y(10)^(1-params(3));
T187 = y(22)^(-1);
T234 = exp(params(23))^(1-params(20));
T245 = exp(params(42))^(1-params(22));
lhs =y(3);
rhs =T10*T16;
residual(1)= lhs-rhs;
lhs =y(1);
rhs =T16-y(3)*params(7)*params(2);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(9)*y(10)/((1-params(3))*y(25));
residual(3)= lhs-rhs;
lhs =y(5)*T38/y(6);
rhs =y(5)*params(2)*T38/y(6)+params(5)/params(10)*(y(4)-(params(5)-1)/params(5));
residual(4)= lhs-rhs;
lhs =y(7)*(y(7)/y(8)-1)/y(8);
rhs =y(7)*params(2)*(y(7)/y(8)-1)/y(8)+T66*T75;
residual(5)= lhs-rhs;
lhs =1;
rhs =T95+y(11);
residual(6)= lhs-rhs;
lhs =y(11);
rhs =exp(x(4))*T89*exp(y(20))*y(13)*params(2)*params(9)-y(13)*exp(y(20))*params(9)*0.5*T103;
residual(7)= lhs-rhs;
lhs =y(17);
rhs =params(43)/params(44)*(exp(params(44)*(y(16)-1))-1);
residual(8)= lhs-rhs;
lhs =y(13);
rhs =T10*params(2)*(y(13)*(1-params(4))+y(16)*y(15)-y(17));
residual(9)= lhs-rhs;
lhs =y(15);
rhs =params(43)*exp(params(44)*(y(16)-1));
residual(10)= lhs-rhs;
lhs =y(15);
rhs =T133/(y(16)*y(18));
residual(11)= lhs-rhs;
lhs =y(9);
rhs =y(9)*y(7)/(y(22)*y(5));
residual(12)= lhs-rhs;
lhs =y(6);
rhs =T144*y(5)^(1-params(16));
residual(13)= lhs-rhs;
lhs =y(8);
rhs =T152*(y(5)*exp(x(4)))^(1-params(18));
residual(14)= lhs-rhs;
lhs =y(1);
rhs =params(2)*exp(y(21))*y(19)*y(1)/(y(22)*y(5));
residual(15)= lhs-rhs;
lhs =y(25);
rhs =T167*(y(2)+y(12)+y(17)*y(18)/y(22));
residual(16)= lhs-rhs;
lhs =y(14);
rhs =1/y(23)-params(10)*0.5*T38^2;
residual(17)= lhs-rhs;
lhs =y(25);
rhs =T182*T183;
residual(18)= lhs-rhs;
lhs =y(18);
rhs =(1-params(4))*y(18)*T187+y(12)*exp(y(20))*(1-T103*params(9)*0.5);
residual(19)= lhs-rhs;
lhs =log(y(19)/params(30));
rhs =log(y(19)/params(30))*params(12)+(1-params(12))*(params(13)*log(y(5)/params(24))+params(14)*y(24)+params(15)*log(y(25)*exp(x(4))/y(25)))+x(5);
residual(20)= lhs-rhs;
lhs =y(24);
rhs =params(3)*log(y(16))+(1-params(3))*log(y(10)/params(25));
residual(21)= lhs-rhs;
lhs =exp(y(23));
rhs =T234*exp(y(23))^params(20)*exp(x(1));
residual(22)= lhs-rhs;
lhs =exp(y(20));
rhs =T245*exp(y(20))^params(22)*exp(x(2));
residual(23)= lhs-rhs;
lhs =y(21);
rhs =y(21)*params(21)+x(3);
residual(24)= lhs-rhs;
lhs =y(22);
rhs =exp(x(4))*params(46);
residual(25)= lhs-rhs;
lhs =y(26);
rhs =y(25);
residual(26)= lhs-rhs;
lhs =y(27);
rhs =y(2);
residual(27)= lhs-rhs;
lhs =y(28);
rhs =y(12);
residual(28)= lhs-rhs;
lhs =y(29);
rhs =x(4);
residual(29)= lhs-rhs;
lhs =y(30);
rhs =exp(x(4))*T89*exp(y(20))*y(13)*params(2)*params(9);
residual(30)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(30, 30);

  %
  % Jacobian matrix
  %

T291 = getPowerDeriv(y(2)-y(2)*params(7)/y(22),(-1),1);
T437 = getPowerDeriv(y(16)*y(18)/y(22),params(3),1);
T463 = 1/params(30)/(y(19)/params(30));
  g1(1,2)=(-(T10*(1-params(7)/y(22))*T291));
  g1(1,3)=1;
  g1(1,22)=(-(T16*(-1)/(y(22)*y(22))+T10*T291*(-((-(y(2)*params(7)))/(y(22)*y(22))))));
  g1(2,1)=1;
  g1(2,2)=(-((1-params(7)/y(22))*T291));
  g1(2,3)=params(7)*params(2);
  g1(2,22)=(-(T291*(-((-(y(2)*params(7)))/(y(22)*y(22))))));
  g1(3,4)=1;
  g1(3,9)=(-(y(10)/((1-params(3))*y(25))));
  g1(3,10)=(-(y(9)/((1-params(3))*y(25))));
  g1(3,25)=(-((-(y(9)*y(10)*(1-params(3))))/((1-params(3))*y(25)*(1-params(3))*y(25))));
  g1(4,4)=(-(params(5)/params(10)));
  g1(4,5)=(T38+y(5)*1/y(6))/y(6)-(params(2)*T38+y(5)*params(2)*1/y(6))/y(6);
  g1(4,6)=(y(6)*y(5)*(-y(5))/(y(6)*y(6))-y(5)*T38)/(y(6)*y(6))-(y(6)*y(5)*params(2)*(-y(5))/(y(6)*y(6))-y(5)*params(2)*T38)/(y(6)*y(6));
  g1(5,1)=(-(T75*y(10)*params(6)/params(11)+T66*(-T70)/(y(1)*y(1))));
  g1(5,7)=(y(7)/y(8)-1+y(7)*1/y(8))/y(8)-(params(2)*(y(7)/y(8)-1)+y(7)*params(2)*1/y(8))/y(8);
  g1(5,8)=(y(8)*y(7)*(-y(7))/(y(8)*y(8))-y(7)*(y(7)/y(8)-1))/(y(8)*y(8))-(y(8)*y(7)*params(2)*(-y(7))/(y(8)*y(8))-y(7)*params(2)*(y(7)/y(8)-1))/(y(8)*y(8));
  g1(5,9)=(-(T66*(-((params(6)-1)/params(6)))));
  g1(5,10)=(-(T75*y(1)*params(6)/params(11)+T66*params(29)*getPowerDeriv(y(10),params(8),1)/y(1)));
  g1(6,11)=(-1);
  g1(6,13)=(-(exp(y(20))*(1-y(12)*exp(x(4))*params(9)*T89/y(12))));
  g1(6,20)=(-T95);
  g1(7,11)=1;
  g1(7,13)=(-(exp(x(4))*T89*exp(y(20))*params(2)*params(9)-T103*0.5*exp(y(20))*params(9)));
  g1(7,20)=(-(exp(x(4))*T89*exp(y(20))*y(13)*params(2)*params(9)-y(13)*exp(y(20))*params(9)*0.5*T103));
  g1(8,16)=(-(params(43)/params(44)*params(44)*exp(params(44)*(y(16)-1))));
  g1(8,17)=1;
  g1(9,13)=1-T10*params(2)*(1-params(4));
  g1(9,15)=(-(y(16)*T10*params(2)));
  g1(9,16)=(-(T10*params(2)*y(15)));
  g1(9,17)=T10*params(2);
  g1(9,22)=(-((y(13)*(1-params(4))+y(16)*y(15)-y(17))*params(2)*(-1)/(y(22)*y(22))));
  g1(10,15)=1;
  g1(10,16)=(-(params(43)*params(44)*exp(params(44)*(y(16)-1))));
  g1(11,9)=(-(y(10)*y(22)*params(3)/(1-params(3))/(y(16)*y(18))));
  g1(11,10)=(-(y(9)*y(22)*params(3)/(1-params(3))/(y(16)*y(18))));
  g1(11,15)=1;
  g1(11,16)=(-((-(T133*y(18)))/(y(16)*y(18)*y(16)*y(18))));
  g1(11,18)=(-((-(y(16)*T133))/(y(16)*y(18)*y(16)*y(18))));
  g1(11,22)=(-(y(10)*y(9)*params(3)/(1-params(3))/(y(16)*y(18))));
  g1(12,5)=(-((-(y(22)*y(9)*y(7)))/(y(22)*y(5)*y(22)*y(5))));
  g1(12,7)=(-(y(9)/(y(22)*y(5))));
  g1(12,9)=1-y(7)/(y(22)*y(5));
  g1(12,22)=(-((-(y(5)*y(9)*y(7)))/(y(22)*y(5)*y(22)*y(5))));
  g1(13,5)=(-(T144*getPowerDeriv(y(5),1-params(16),1)));
  g1(13,6)=1;
  g1(14,5)=(-(T152*exp(x(4))*getPowerDeriv(y(5)*exp(x(4)),1-params(18),1)));
  g1(14,8)=1;
  g1(15,1)=1-params(2)*exp(y(21))*y(19)*1/(y(22)*y(5));
  g1(15,5)=(-(params(2)*exp(y(21))*y(19)*(-(y(22)*y(1)))/(y(22)*y(5)*y(22)*y(5))));
  g1(15,19)=(-(params(2)*exp(y(21))*y(1)/(y(22)*y(5))));
  g1(15,21)=(-(params(2)*exp(y(21))*y(19)*y(1)/(y(22)*y(5))));
  g1(15,22)=(-(params(2)*exp(y(21))*y(19)*(-(y(1)*y(5)))/(y(22)*y(5)*y(22)*y(5))));
  g1(16,2)=(-T167);
  g1(16,12)=(-T167);
  g1(16,14)=(-((y(2)+y(12)+y(17)*y(18)/y(22))*getPowerDeriv(y(14),(-1),1)));
  g1(16,17)=(-(T167*y(18)/y(22)));
  g1(16,18)=(-(T167*y(17)/y(22)));
  g1(16,22)=(-(T167*(-(y(17)*y(18)))/(y(22)*y(22))));
  g1(16,25)=1;
  g1(17,5)=params(10)*0.5*1/y(6)*2*T38;
  g1(17,6)=params(10)*0.5*2*T38*(-y(5))/(y(6)*y(6));
  g1(17,14)=1;
  g1(17,23)=(-((-1)/(y(23)*y(23))));
  g1(18,10)=(-(T182*getPowerDeriv(y(10),1-params(3),1)));
  g1(18,16)=(-(T183*y(18)/y(22)*T437));
  g1(18,18)=(-(T183*T437*y(16)/y(22)));
  g1(18,22)=(-(T183*T437*(-(y(16)*y(18)))/(y(22)*y(22))));
  g1(18,25)=1;
  g1(19,12)=(-(exp(y(20))*(1-T103*params(9)*0.5)));
  g1(19,18)=1-(1-params(4))*T187;
  g1(19,20)=(-(y(12)*exp(y(20))*(1-T103*params(9)*0.5)));
  g1(19,22)=(-((1-params(4))*y(18)*getPowerDeriv(y(22),(-1),1)));
  g1(20,5)=(-((1-params(12))*params(13)*1/params(24)/(y(5)/params(24))));
  g1(20,19)=T463-params(12)*T463;
  g1(20,24)=(-((1-params(12))*params(14)));
  g1(21,10)=(-((1-params(3))*1/params(25)/(y(10)/params(25))));
  g1(21,16)=(-(params(3)*1/y(16)));
  g1(21,24)=1;
  g1(22,23)=exp(y(23))-exp(x(1))*T234*exp(y(23))*getPowerDeriv(exp(y(23)),params(20),1);
  g1(23,20)=exp(y(20))-exp(x(2))*T245*exp(y(20))*getPowerDeriv(exp(y(20)),params(22),1);
  g1(24,21)=1-params(21);
  g1(25,22)=1;
  g1(26,25)=(-1);
  g1(26,26)=1;
  g1(27,2)=(-1);
  g1(27,27)=1;
  g1(28,12)=(-1);
  g1(28,28)=1;
  g1(29,29)=1;
  g1(30,13)=(-(exp(x(4))*T89*exp(y(20))*params(2)*params(9)));
  g1(30,20)=(-(exp(x(4))*T89*exp(y(20))*y(13)*params(2)*params(9)));
  g1(30,30)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],30,900);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],30,27000);
end
end
end
end

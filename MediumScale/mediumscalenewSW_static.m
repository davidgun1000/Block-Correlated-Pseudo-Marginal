function [residual, g1, g2, g3] = mediumscalenewSW_static(y, x, params)
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

residual = zeros( 51, 1);

%
% Model equations
%

T10 = 1/y(22);
T16 = (y(2)-y(2)*params(7)/y(22))^(-1);
T39 = y(5)/y(6)-1;
T72 = params(29)*y(10)^params(8);
T77 = T72/y(1)-y(9)*(params(6)-1)/params(6);
T91 = y(12)*exp(x(4))/y(12)-1;
T97 = y(13)*exp(y(20))*(1-y(12)*exp(x(4))*params(9)*T91/y(12));
T105 = T91^2;
T112 = params(43)/params(44);
T133 = y(22)*params(3)/(1-params(3));
T135 = y(10)*y(9)*T133;
T146 = params(24)^params(16);
T154 = params(46)*params(24)^params(18);
T169 = y(14)^(-1);
T184 = (y(16)*y(18)/y(22))^params(3);
T185 = y(10)^(1-params(3));
T189 = y(22)^(-1);
T239 = exp(params(23))^(1-params(20));
T250 = exp(params(42))^(1-params(22));
T277 = (y(31)-params(7)*y(31)/y(22))^(-1);
T322 = params(29)*y(35)^params(8);
T326 = T322/y(32)-(params(6)-1)*y(34)/params(6);
T335 = params(50)^params(16);
T340 = params(46)*params(50)^params(18);
T358 = exp(x(4))*y(38)/y(38)-1;
T364 = exp(y(20))*y(37)*(1-y(38)*exp(x(4))*params(9)*T358/y(38));
T371 = T358^2;
T396 = y(44)^(-1);
T405 = (y(41)*y(45)/y(22))^params(3);
T406 = y(35)^(1-params(3));
lhs =y(3);
rhs =T10*T16;
residual(1)= lhs-rhs;
lhs =y(1);
rhs =T16-y(3)*params(7)*params(2);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(9)*y(10)/((1-params(3))*y(25));
residual(3)= lhs-rhs;
lhs =y(5)*params(10)*T39/y(6);
rhs =y(5)*T39*params(2)*params(10)/y(6)+params(5)*(y(4)-(params(5)-1)/params(5));
residual(4)= lhs-rhs;
lhs =y(7)*params(11)*(y(7)/y(8)-1)/y(8);
rhs =y(7)*(y(7)/y(8)-1)*params(2)*params(11)/y(8)+y(1)*y(10)*params(6)*T77;
residual(5)= lhs-rhs;
lhs =1;
rhs =T97+y(11);
residual(6)= lhs-rhs;
lhs =y(11);
rhs =exp(x(4))*T91*exp(y(20))*y(13)*params(2)*params(9)-y(13)*exp(y(20))*params(9)*0.5*T105;
residual(7)= lhs-rhs;
lhs =y(17);
rhs =T112*(exp(params(44)*(y(16)-1))-1);
residual(8)= lhs-rhs;
lhs =y(13);
rhs =T10*params(2)*(y(13)*(1-params(4))+y(16)*y(15)-y(17));
residual(9)= lhs-rhs;
lhs =y(15);
rhs =params(43)*exp(params(44)*(y(16)-1));
residual(10)= lhs-rhs;
lhs =y(15);
rhs =T135/(y(16)*y(18));
residual(11)= lhs-rhs;
lhs =y(9);
rhs =y(9)*y(7)/(y(22)*y(5));
residual(12)= lhs-rhs;
lhs =y(6);
rhs =T146*y(5)^(1-params(16));
residual(13)= lhs-rhs;
lhs =y(8);
rhs =T154*(y(5)*exp(x(4)))^(1-params(18));
residual(14)= lhs-rhs;
lhs =y(1);
rhs =params(2)*exp(y(21))*y(19)*y(1)/(y(22)*y(5));
residual(15)= lhs-rhs;
lhs =y(25);
rhs =T169*(y(2)+y(12)+y(17)*y(18)/y(22));
residual(16)= lhs-rhs;
lhs =y(14);
rhs =1/y(23)-params(10)*0.5*T39^2;
residual(17)= lhs-rhs;
lhs =y(25);
rhs =T184*T185;
residual(18)= lhs-rhs;
lhs =y(18);
rhs =(1-params(4))*y(18)*T189+y(12)*exp(y(20))*(1-T105*params(9)*0.5);
residual(19)= lhs-rhs;
lhs =log(y(19)/params(30));
rhs =log(y(19)/params(30))*params(12)+(1-params(12))*(params(13)*log(y(5)/params(24))+params(14)*log(y(25)/y(36))+params(15)*log(y(25)*exp(x(4))/y(25)))+x(5);
residual(20)= lhs-rhs;
lhs =y(24);
rhs =params(3)*log(y(16))+(1-params(3))*log(y(10)/params(25));
residual(21)= lhs-rhs;
lhs =exp(y(23));
rhs =T239*exp(y(23))^params(20)*exp(x(1));
residual(22)= lhs-rhs;
lhs =exp(y(20));
rhs =T250*exp(y(20))^params(22)*exp(x(2));
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
rhs =T10*T277;
residual(30)= lhs-rhs;
lhs =y(32);
rhs =T277-params(7)*params(2)*y(30);
residual(31)= lhs-rhs;
lhs =y(33);
rhs =y(34)*y(35)/((1-params(3))*y(36));
residual(32)= lhs-rhs;
lhs =y(46)*params(48)*(y(46)/y(47)-1)/y(47);
rhs =y(46)*(y(46)/y(47)-1)*params(2)*params(48)/y(47)+params(5)*(y(33)-(params(5)-1)/params(5));
residual(33)= lhs-rhs;
lhs =y(48)*params(49)*(y(48)/y(49)-1)/y(49);
rhs =y(48)*(y(48)/y(49)-1)*params(2)*params(49)/y(49)+params(6)*y(32)*y(35)*T326;
residual(34)= lhs-rhs;
lhs =y(34);
rhs =y(34)*y(48)/(y(22)*y(46));
residual(35)= lhs-rhs;
lhs =y(47);
rhs =T335*y(46)^(1-params(16));
residual(36)= lhs-rhs;
lhs =y(49);
rhs =T340*(exp(x(4))*y(46))^(1-params(18));
residual(37)= lhs-rhs;
lhs =y(42);
rhs =y(35)*T133*y(34)/(y(41)*y(45));
residual(38)= lhs-rhs;
lhs =1;
rhs =T364+y(39);
residual(39)= lhs-rhs;
lhs =y(39);
rhs =exp(x(4))*T358*exp(y(20))*y(37)*params(2)*params(9)-0.5*params(9)*exp(y(20))*y(37)*T371;
residual(40)= lhs-rhs;
lhs =y(40);
rhs =T112*(exp(params(44)*(y(41)-1))-1);
residual(41)= lhs-rhs;
lhs =y(37);
rhs =T10*params(2)*((1-params(4))*y(37)+y(42)*y(41)-y(40));
residual(42)= lhs-rhs;
lhs =y(42);
rhs =params(43)*exp(params(44)*(y(41)-1));
residual(43)= lhs-rhs;
lhs =y(32);
rhs =params(2)*exp(y(21))*y(43)*y(32)/(y(22)*y(46));
residual(44)= lhs-rhs;
lhs =y(36);
rhs =T396*(y(31)+y(38)+y(45)*y(40)/y(22));
residual(45)= lhs-rhs;
lhs =y(44);
rhs =1/y(23);
residual(46)= lhs-rhs;
lhs =y(36);
rhs =T405*T406;
residual(47)= lhs-rhs;
lhs =y(45);
rhs =T189*(1-params(4))*y(45)+y(38)*exp(y(20))*(1-params(9)*0.5*T371);
residual(48)= lhs-rhs;
lhs =log(y(43)/params(30));
rhs =x(5)+params(12)*log(y(43)/params(30))+(1-params(12))*(params(13)*log(y(46)/params(24))+params(15)*log(exp(x(4))*y(36)/y(36)));
residual(49)= lhs-rhs;
lhs =y(50);
rhs =exp(x(4))*T91*exp(y(20))*y(13)*params(2)*params(9);
residual(50)= lhs-rhs;
lhs =y(51);
rhs =exp(x(4))*T358*exp(y(20))*y(37)*params(2)*params(9);
residual(51)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(51, 51);

  %
  % Jacobian matrix
  %

T466 = getPowerDeriv(y(2)-y(2)*params(7)/y(22),(-1),1);
T616 = getPowerDeriv(y(16)*y(18)/y(22),params(3),1);
T642 = 1/params(30)/(y(19)/params(30));
T697 = getPowerDeriv(y(22),(-1),1);
T703 = getPowerDeriv(y(31)-params(7)*y(31)/y(22),(-1),1);
T732 = getPowerDeriv(y(41)*y(45)/y(22),params(3),1);
T841 = 1/params(30)/(y(43)/params(30));
  g1(1,2)=(-(T10*(1-params(7)/y(22))*T466));
  g1(1,3)=1;
  g1(1,22)=(-(T16*(-1)/(y(22)*y(22))+T10*T466*(-((-(y(2)*params(7)))/(y(22)*y(22))))));
  g1(2,1)=1;
  g1(2,2)=(-((1-params(7)/y(22))*T466));
  g1(2,3)=params(7)*params(2);
  g1(2,22)=(-(T466*(-((-(y(2)*params(7)))/(y(22)*y(22))))));
  g1(3,4)=1;
  g1(3,9)=(-(y(10)/((1-params(3))*y(25))));
  g1(3,10)=(-(y(9)/((1-params(3))*y(25))));
  g1(3,25)=(-((-(y(9)*y(10)*(1-params(3))))/((1-params(3))*y(25)*(1-params(3))*y(25))));
  g1(4,4)=(-params(5));
  g1(4,5)=(params(10)*T39+y(5)*params(10)*1/y(6))/y(6)-(T39*params(2)*params(10)+y(5)*params(2)*params(10)*1/y(6))/y(6);
  g1(4,6)=(y(6)*y(5)*params(10)*(-y(5))/(y(6)*y(6))-y(5)*params(10)*T39)/(y(6)*y(6))-(y(6)*y(5)*params(2)*params(10)*(-y(5))/(y(6)*y(6))-y(5)*T39*params(2)*params(10))/(y(6)*y(6));
  g1(5,1)=(-(T77*y(10)*params(6)+y(1)*y(10)*params(6)*(-T72)/(y(1)*y(1))));
  g1(5,7)=(params(11)*(y(7)/y(8)-1)+y(7)*params(11)*1/y(8))/y(8)-((y(7)/y(8)-1)*params(2)*params(11)+y(7)*params(2)*params(11)*1/y(8))/y(8);
  g1(5,8)=(y(8)*y(7)*params(11)*(-y(7))/(y(8)*y(8))-y(7)*params(11)*(y(7)/y(8)-1))/(y(8)*y(8))-(y(8)*y(7)*params(2)*params(11)*(-y(7))/(y(8)*y(8))-y(7)*(y(7)/y(8)-1)*params(2)*params(11))/(y(8)*y(8));
  g1(5,9)=(-(y(1)*y(10)*params(6)*(-((params(6)-1)/params(6)))));
  g1(5,10)=(-(T77*y(1)*params(6)+y(1)*y(10)*params(6)*params(29)*getPowerDeriv(y(10),params(8),1)/y(1)));
  g1(6,11)=(-1);
  g1(6,13)=(-(exp(y(20))*(1-y(12)*exp(x(4))*params(9)*T91/y(12))));
  g1(6,20)=(-T97);
  g1(7,11)=1;
  g1(7,13)=(-(exp(x(4))*T91*exp(y(20))*params(2)*params(9)-T105*0.5*exp(y(20))*params(9)));
  g1(7,20)=(-(exp(x(4))*T91*exp(y(20))*y(13)*params(2)*params(9)-y(13)*exp(y(20))*params(9)*0.5*T105));
  g1(8,16)=(-(T112*params(44)*exp(params(44)*(y(16)-1))));
  g1(8,17)=1;
  g1(9,13)=1-T10*params(2)*(1-params(4));
  g1(9,15)=(-(y(16)*T10*params(2)));
  g1(9,16)=(-(T10*params(2)*y(15)));
  g1(9,17)=T10*params(2);
  g1(9,22)=(-((y(13)*(1-params(4))+y(16)*y(15)-y(17))*params(2)*(-1)/(y(22)*y(22))));
  g1(10,15)=1;
  g1(10,16)=(-(params(43)*params(44)*exp(params(44)*(y(16)-1))));
  g1(11,9)=(-(y(10)*T133/(y(16)*y(18))));
  g1(11,10)=(-(y(9)*T133/(y(16)*y(18))));
  g1(11,15)=1;
  g1(11,16)=(-((-(T135*y(18)))/(y(16)*y(18)*y(16)*y(18))));
  g1(11,18)=(-((-(y(16)*T135))/(y(16)*y(18)*y(16)*y(18))));
  g1(11,22)=(-(y(10)*y(9)*params(3)/(1-params(3))/(y(16)*y(18))));
  g1(12,5)=(-((-(y(22)*y(9)*y(7)))/(y(22)*y(5)*y(22)*y(5))));
  g1(12,7)=(-(y(9)/(y(22)*y(5))));
  g1(12,9)=1-y(7)/(y(22)*y(5));
  g1(12,22)=(-((-(y(5)*y(9)*y(7)))/(y(22)*y(5)*y(22)*y(5))));
  g1(13,5)=(-(T146*getPowerDeriv(y(5),1-params(16),1)));
  g1(13,6)=1;
  g1(14,5)=(-(T154*exp(x(4))*getPowerDeriv(y(5)*exp(x(4)),1-params(18),1)));
  g1(14,8)=1;
  g1(15,1)=1-params(2)*exp(y(21))*y(19)*1/(y(22)*y(5));
  g1(15,5)=(-(params(2)*exp(y(21))*y(19)*(-(y(22)*y(1)))/(y(22)*y(5)*y(22)*y(5))));
  g1(15,19)=(-(params(2)*exp(y(21))*y(1)/(y(22)*y(5))));
  g1(15,21)=(-(params(2)*exp(y(21))*y(19)*y(1)/(y(22)*y(5))));
  g1(15,22)=(-(params(2)*exp(y(21))*y(19)*(-(y(1)*y(5)))/(y(22)*y(5)*y(22)*y(5))));
  g1(16,2)=(-T169);
  g1(16,12)=(-T169);
  g1(16,14)=(-((y(2)+y(12)+y(17)*y(18)/y(22))*getPowerDeriv(y(14),(-1),1)));
  g1(16,17)=(-(T169*y(18)/y(22)));
  g1(16,18)=(-(T169*y(17)/y(22)));
  g1(16,22)=(-(T169*(-(y(17)*y(18)))/(y(22)*y(22))));
  g1(16,25)=1;
  g1(17,5)=params(10)*0.5*1/y(6)*2*T39;
  g1(17,6)=params(10)*0.5*2*T39*(-y(5))/(y(6)*y(6));
  g1(17,14)=1;
  g1(17,23)=(-((-1)/(y(23)*y(23))));
  g1(18,10)=(-(T184*getPowerDeriv(y(10),1-params(3),1)));
  g1(18,16)=(-(T185*y(18)/y(22)*T616));
  g1(18,18)=(-(T185*T616*y(16)/y(22)));
  g1(18,22)=(-(T185*T616*(-(y(16)*y(18)))/(y(22)*y(22))));
  g1(18,25)=1;
  g1(19,12)=(-(exp(y(20))*(1-T105*params(9)*0.5)));
  g1(19,18)=1-(1-params(4))*T189;
  g1(19,20)=(-(y(12)*exp(y(20))*(1-T105*params(9)*0.5)));
  g1(19,22)=(-((1-params(4))*y(18)*T697));
  g1(20,5)=(-((1-params(12))*params(13)*1/params(24)/(y(5)/params(24))));
  g1(20,19)=T642-params(12)*T642;
  g1(20,25)=(-((1-params(12))*params(14)*1/y(36)/(y(25)/y(36))));
  g1(20,36)=(-((1-params(12))*params(14)*(-y(25))/(y(36)*y(36))/(y(25)/y(36))));
  g1(21,10)=(-((1-params(3))*1/params(25)/(y(10)/params(25))));
  g1(21,16)=(-(params(3)*1/y(16)));
  g1(21,24)=1;
  g1(22,23)=exp(y(23))-exp(x(1))*T239*exp(y(23))*getPowerDeriv(exp(y(23)),params(20),1);
  g1(23,20)=exp(y(20))-exp(x(2))*T250*exp(y(20))*getPowerDeriv(exp(y(20)),params(22),1);
  g1(24,21)=1-params(21);
  g1(25,22)=1;
  g1(26,25)=(-1);
  g1(26,26)=1;
  g1(27,2)=(-1);
  g1(27,27)=1;
  g1(28,12)=(-1);
  g1(28,28)=1;
  g1(29,29)=1;
  g1(30,22)=(-(T277*(-1)/(y(22)*y(22))+T10*(-((-(params(7)*y(31)))/(y(22)*y(22))))*T703));
  g1(30,30)=1;
  g1(30,31)=(-(T10*(1-params(7)/y(22))*T703));
  g1(31,22)=(-((-((-(params(7)*y(31)))/(y(22)*y(22))))*T703));
  g1(31,30)=params(7)*params(2);
  g1(31,31)=(-((1-params(7)/y(22))*T703));
  g1(31,32)=1;
  g1(32,33)=1;
  g1(32,34)=(-(y(35)/((1-params(3))*y(36))));
  g1(32,35)=(-(y(34)/((1-params(3))*y(36))));
  g1(32,36)=(-((-((1-params(3))*y(34)*y(35)))/((1-params(3))*y(36)*(1-params(3))*y(36))));
  g1(33,33)=(-params(5));
  g1(33,46)=(params(48)*(y(46)/y(47)-1)+y(46)*params(48)*1/y(47))/y(47)-((y(46)/y(47)-1)*params(2)*params(48)+y(46)*params(2)*params(48)*1/y(47))/y(47);
  g1(33,47)=(y(47)*y(46)*params(48)*(-y(46))/(y(47)*y(47))-y(46)*params(48)*(y(46)/y(47)-1))/(y(47)*y(47))-(y(47)*y(46)*params(2)*params(48)*(-y(46))/(y(47)*y(47))-y(46)*(y(46)/y(47)-1)*params(2)*params(48))/(y(47)*y(47));
  g1(34,32)=(-(T326*params(6)*y(35)+params(6)*y(32)*y(35)*(-T322)/(y(32)*y(32))));
  g1(34,34)=(-(params(6)*y(32)*y(35)*(-((params(6)-1)/params(6)))));
  g1(34,35)=(-(T326*params(6)*y(32)+params(6)*y(32)*y(35)*params(29)*getPowerDeriv(y(35),params(8),1)/y(32)));
  g1(34,48)=(params(49)*(y(48)/y(49)-1)+y(48)*params(49)*1/y(49))/y(49)-((y(48)/y(49)-1)*params(2)*params(49)+y(48)*params(2)*params(49)*1/y(49))/y(49);
  g1(34,49)=(y(49)*y(48)*params(49)*(-y(48))/(y(49)*y(49))-y(48)*params(49)*(y(48)/y(49)-1))/(y(49)*y(49))-(y(49)*y(48)*params(2)*params(49)*(-y(48))/(y(49)*y(49))-y(48)*(y(48)/y(49)-1)*params(2)*params(49))/(y(49)*y(49));
  g1(35,22)=(-((-(y(46)*y(34)*y(48)))/(y(22)*y(46)*y(22)*y(46))));
  g1(35,34)=1-y(48)/(y(22)*y(46));
  g1(35,46)=(-((-(y(22)*y(34)*y(48)))/(y(22)*y(46)*y(22)*y(46))));
  g1(35,48)=(-(y(34)/(y(22)*y(46))));
  g1(36,46)=(-(T335*getPowerDeriv(y(46),1-params(16),1)));
  g1(36,47)=1;
  g1(37,46)=(-(T340*exp(x(4))*getPowerDeriv(exp(x(4))*y(46),1-params(18),1)));
  g1(37,49)=1;
  g1(38,22)=(-(y(35)*params(3)/(1-params(3))*y(34)/(y(41)*y(45))));
  g1(38,34)=(-(T133*y(35)/(y(41)*y(45))));
  g1(38,35)=(-(T133*y(34)/(y(41)*y(45))));
  g1(38,41)=(-((-(y(35)*T133*y(34)*y(45)))/(y(41)*y(45)*y(41)*y(45))));
  g1(38,42)=1;
  g1(38,45)=(-((-(y(35)*T133*y(34)*y(41)))/(y(41)*y(45)*y(41)*y(45))));
  g1(39,20)=(-T364);
  g1(39,37)=(-(exp(y(20))*(1-y(38)*exp(x(4))*params(9)*T358/y(38))));
  g1(39,39)=(-1);
  g1(40,20)=(-(exp(x(4))*T358*exp(y(20))*y(37)*params(2)*params(9)-0.5*params(9)*exp(y(20))*y(37)*T371));
  g1(40,37)=(-(exp(x(4))*T358*exp(y(20))*params(2)*params(9)-T371*0.5*exp(y(20))*params(9)));
  g1(40,39)=1;
  g1(41,40)=1;
  g1(41,41)=(-(T112*params(44)*exp(params(44)*(y(41)-1))));
  g1(42,22)=(-(((1-params(4))*y(37)+y(42)*y(41)-y(40))*params(2)*(-1)/(y(22)*y(22))));
  g1(42,37)=1-T10*params(2)*(1-params(4));
  g1(42,40)=T10*params(2);
  g1(42,41)=(-(T10*params(2)*y(42)));
  g1(42,42)=(-(T10*params(2)*y(41)));
  g1(43,41)=(-(params(43)*params(44)*exp(params(44)*(y(41)-1))));
  g1(43,42)=1;
  g1(44,21)=(-(params(2)*exp(y(21))*y(43)*y(32)/(y(22)*y(46))));
  g1(44,22)=(-(params(2)*exp(y(21))*y(43)*(-(y(32)*y(46)))/(y(22)*y(46)*y(22)*y(46))));
  g1(44,32)=1-params(2)*exp(y(21))*y(43)*1/(y(22)*y(46));
  g1(44,43)=(-(params(2)*exp(y(21))*y(32)/(y(22)*y(46))));
  g1(44,46)=(-(params(2)*exp(y(21))*y(43)*(-(y(22)*y(32)))/(y(22)*y(46)*y(22)*y(46))));
  g1(45,22)=(-(T396*(-(y(45)*y(40)))/(y(22)*y(22))));
  g1(45,31)=(-T396);
  g1(45,36)=1;
  g1(45,38)=(-T396);
  g1(45,40)=(-(T396*y(45)/y(22)));
  g1(45,44)=(-((y(31)+y(38)+y(45)*y(40)/y(22))*getPowerDeriv(y(44),(-1),1)));
  g1(45,45)=(-(T396*y(40)/y(22)));
  g1(46,23)=(-((-1)/(y(23)*y(23))));
  g1(46,44)=1;
  g1(47,22)=(-(T406*(-(y(41)*y(45)))/(y(22)*y(22))*T732));
  g1(47,35)=(-(T405*getPowerDeriv(y(35),1-params(3),1)));
  g1(47,36)=1;
  g1(47,41)=(-(T406*T732*y(45)/y(22)));
  g1(47,45)=(-(T406*T732*y(41)/y(22)));
  g1(48,20)=(-(y(38)*exp(y(20))*(1-params(9)*0.5*T371)));
  g1(48,22)=(-((1-params(4))*y(45)*T697));
  g1(48,38)=(-(exp(y(20))*(1-params(9)*0.5*T371)));
  g1(48,45)=1-(1-params(4))*T189;
  g1(49,43)=T841-params(12)*T841;
  g1(49,46)=(-((1-params(12))*params(13)*1/params(24)/(y(46)/params(24))));
  g1(50,13)=(-(exp(x(4))*T91*exp(y(20))*params(2)*params(9)));
  g1(50,20)=(-(exp(x(4))*T91*exp(y(20))*y(13)*params(2)*params(9)));
  g1(50,50)=1;
  g1(51,20)=(-(exp(x(4))*T358*exp(y(20))*y(37)*params(2)*params(9)));
  g1(51,37)=(-(exp(x(4))*T358*exp(y(20))*params(2)*params(9)));
  g1(51,51)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],51,2601);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],51,132651);
end
end
end
end

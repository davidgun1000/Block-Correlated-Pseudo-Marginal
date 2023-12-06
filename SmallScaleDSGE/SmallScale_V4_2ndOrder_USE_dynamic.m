function [residual, g1, g2, g3] = SmallScale_V4_2ndOrder_USE_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(8, 1);
T13 = exp(y(14))/exp(y(7));
T16 = T13^(-params(5));
T25 = exp(y(15))*exp(y(16));
T26 = params(1)*T16*exp(y(8))/T25;
T39 = exp(y(10))*(1-0.5/params(2));
T42 = T39+exp(params(4))*0.5/params(2);
T47 = exp(y(15))*params(1)*params(3)*(exp(y(15))-exp(params(4)));
T53 = exp(y(13))/exp(y(5));
T54 = T16*T47*T53;
T66 = exp(y(9))^(1-params(9));
T69 = exp(y(2))^params(9);
T77 = exp(params(4))*params(13)/params(1);
T78 = exp(y(10))/exp(params(4));
T81 = T77*T78^params(6);
T87 = (1-params(2))^(1/params(5))*exp(y(11));
T88 = exp(y(5))/T87;
T90 = T88^params(7);
T93 = exp(y(7))/exp(y(5));
T97 = params(3)*0.5*exp(params(4))^2;
lhs =T26;
rhs =1;
residual(1)= lhs-rhs;
lhs =params(3)*(exp(y(10))-exp(params(4)))*T42-T54+1/params(2)*(1-exp(y(7))^params(5));
rhs =1;
residual(2)= lhs-rhs;
lhs =exp(y(8));
rhs =T66*T69*exp(x(it_, 1));
residual(3)= lhs-rhs;
lhs =exp(y(9));
rhs =T81*T90;
residual(4)= lhs-rhs;
lhs =T93;
rhs =1/exp(y(11))-T97*(T78-1)^2;
residual(5)= lhs-rhs;
lhs =y(11);
rhs =(1-params(8))*params(10)+params(8)*y(3)+x(it_, 2);
residual(6)= lhs-rhs;
lhs =y(12);
rhs =(1-params(17))*params(14)+params(17)*y(4)+x(it_, 3);
residual(7)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(1));
residual(8)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(8, 19);

  %
  % Jacobian matrix
  %

T133 = (-(exp(y(13))*exp(y(5))))/(exp(y(5))*exp(y(5)));
T135 = (-(T16*T47*T133));
T136 = getPowerDeriv(T88,params(7),1);
T142 = (-(exp(y(7))*exp(y(5))))/(exp(y(5))*exp(y(5)));
T147 = (-(exp(y(14))*exp(y(7))))/(exp(y(7))*exp(y(7)));
T148 = getPowerDeriv(T13,(-params(5)),1);
T149 = T147*T148;
T157 = exp(y(7))*getPowerDeriv(exp(y(7)),params(5),1);
T169 = exp(y(2))*getPowerDeriv(exp(y(2)),params(9),1);
T174 = exp(y(9))*getPowerDeriv(exp(y(9)),1-params(9),1);
T183 = T78*getPowerDeriv(T78,params(6),1);
T184 = T77*T183;
T194 = (-(params(1)*T16*exp(y(8))*T25))/(T25*T25);
T205 = (-(exp(y(5))*T87))/(T87*T87);
T206 = T136*T205;
  g1(1,7)=exp(y(8))*params(1)*T149/T25;
  g1(1,14)=exp(y(8))*params(1)*T13*T148/T25;
  g1(1,8)=T26;
  g1(1,15)=T194;
  g1(1,16)=T194;
  g1(2,5)=T135;
  g1(2,13)=(-T54);
  g1(2,7)=(-(T53*T47*T149))+1/params(2)*(-T157);
  g1(2,14)=(-(T53*T47*T13*T148));
  g1(2,10)=T42*params(3)*exp(y(10))+params(3)*(exp(y(10))-exp(params(4)))*T39;
  g1(2,15)=(-(T53*T16*(T47+exp(y(15))*exp(y(15))*params(1)*params(3))));
  g1(3,2)=(-(exp(x(it_, 1))*T66*T169));
  g1(3,8)=exp(y(8));
  g1(3,9)=(-(exp(x(it_, 1))*T69*T174));
  g1(3,17)=(-(T66*T69*exp(x(it_, 1))));
  g1(4,5)=(-(T81*T88*T136));
  g1(4,9)=exp(y(9));
  g1(4,10)=(-(T90*T184));
  g1(4,11)=(-(T81*T206));
  g1(5,5)=T142;
  g1(5,7)=T93;
  g1(5,10)=T97*T78*2*(T78-1);
  g1(5,11)=(-((-exp(y(11)))/(exp(y(11))*exp(y(11)))));
  g1(6,3)=(-params(8));
  g1(6,11)=1;
  g1(6,18)=(-1);
  g1(7,4)=(-params(17));
  g1(7,12)=1;
  g1(7,19)=(-1);
  g1(8,1)=(-exp(y(1)));
  g1(8,6)=exp(y(6));

if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(79,3);
T221 = getPowerDeriv(T13,(-params(5)),2);
T222 = T147*T221;
T225 = T148*((-(exp(y(14))*exp(y(7))))*exp(y(7))*exp(y(7))-(-(exp(y(14))*exp(y(7))))*(exp(y(7))*exp(y(7))+exp(y(7))*exp(y(7))))/(exp(y(7))*exp(y(7))*exp(y(7))*exp(y(7)))+T147*T222;
T251 = ((-(params(1)*T16*exp(y(8))*T25))*T25*T25-(-(params(1)*T16*exp(y(8))*T25))*(T25*T25+T25*T25))/(T25*T25*T25*T25);
T314 = getPowerDeriv(T88,params(7),2);
T315 = T88*T314;
  v2(1,1)=1;
  v2(1,2)=121;
  v2(1,3)=exp(y(8))*params(1)*T225/T25;
  v2(2,1)=1;
  v2(2,2)=254;
  v2(2,3)=exp(y(8))*params(1)*(T149+T13*T222)/T25;
  v2(3,1)=1;
  v2(3,2)=128;
  v2(3,3)=  v2(2,3);
  v2(4,1)=1;
  v2(4,2)=261;
  v2(4,3)=exp(y(8))*params(1)*(T13*T148+T13*T13*T221)/T25;
  v2(5,1)=1;
  v2(5,2)=140;
  v2(5,3)=exp(y(8))*params(1)*T149/T25;
  v2(6,1)=1;
  v2(6,2)=122;
  v2(6,3)=  v2(5,3);
  v2(7,1)=1;
  v2(7,2)=147;
  v2(7,3)=exp(y(8))*params(1)*T13*T148/T25;
  v2(8,1)=1;
  v2(8,2)=255;
  v2(8,3)=  v2(7,3);
  v2(9,1)=1;
  v2(9,2)=141;
  v2(9,3)=T26;
  v2(10,1)=1;
  v2(10,2)=273;
  v2(10,3)=(-(T25*exp(y(8))*params(1)*T149))/(T25*T25);
  v2(11,1)=1;
  v2(11,2)=129;
  v2(11,3)=  v2(10,3);
  v2(12,1)=1;
  v2(12,2)=280;
  v2(12,3)=(-(T25*exp(y(8))*params(1)*T13*T148))/(T25*T25);
  v2(13,1)=1;
  v2(13,2)=262;
  v2(13,3)=  v2(12,3);
  v2(14,1)=1;
  v2(14,2)=274;
  v2(14,3)=T194;
  v2(15,1)=1;
  v2(15,2)=148;
  v2(15,3)=  v2(14,3);
  v2(16,1)=1;
  v2(16,2)=281;
  v2(16,3)=T251;
  v2(17,1)=1;
  v2(17,2)=292;
  v2(17,3)=(-(T25*exp(y(8))*params(1)*T149))/(T25*T25);
  v2(18,1)=1;
  v2(18,2)=130;
  v2(18,3)=  v2(17,3);
  v2(19,1)=1;
  v2(19,2)=299;
  v2(19,3)=(-(T25*exp(y(8))*params(1)*T13*T148))/(T25*T25);
  v2(20,1)=1;
  v2(20,2)=263;
  v2(20,3)=  v2(19,3);
  v2(21,1)=1;
  v2(21,2)=293;
  v2(21,3)=T194;
  v2(22,1)=1;
  v2(22,2)=149;
  v2(22,3)=  v2(21,3);
  v2(23,1)=1;
  v2(23,2)=300;
  v2(23,3)=T251;
  v2(24,1)=1;
  v2(24,2)=282;
  v2(24,3)=  v2(23,3);
  v2(25,1)=1;
  v2(25,2)=301;
  v2(25,3)=T251;
  v2(26,1)=2;
  v2(26,2)=81;
  v2(26,3)=(-(T16*T47*((-(exp(y(13))*exp(y(5))))*exp(y(5))*exp(y(5))-(-(exp(y(13))*exp(y(5))))*(exp(y(5))*exp(y(5))+exp(y(5))*exp(y(5))))/(exp(y(5))*exp(y(5))*exp(y(5))*exp(y(5)))));
  v2(27,1)=2;
  v2(27,2)=233;
  v2(27,3)=T135;
  v2(28,1)=2;
  v2(28,2)=89;
  v2(28,3)=  v2(27,3);
  v2(29,1)=2;
  v2(29,2)=241;
  v2(29,3)=(-T54);
  v2(30,1)=2;
  v2(30,2)=119;
  v2(30,3)=(-(T133*T47*T149));
  v2(31,1)=2;
  v2(31,2)=83;
  v2(31,3)=  v2(30,3);
  v2(32,1)=2;
  v2(32,2)=127;
  v2(32,3)=(-(T53*T47*T149));
  v2(33,1)=2;
  v2(33,2)=235;
  v2(33,3)=  v2(32,3);
  v2(34,1)=2;
  v2(34,2)=121;
  v2(34,3)=(-(T53*T47*T225))+1/params(2)*(-(T157+exp(y(7))*exp(y(7))*getPowerDeriv(exp(y(7)),params(5),2)));
  v2(35,1)=2;
  v2(35,2)=252;
  v2(35,3)=(-(T133*T47*T13*T148));
  v2(36,1)=2;
  v2(36,2)=90;
  v2(36,3)=  v2(35,3);
  v2(37,1)=2;
  v2(37,2)=260;
  v2(37,3)=(-(T53*T47*T13*T148));
  v2(38,1)=2;
  v2(38,2)=242;
  v2(38,3)=  v2(37,3);
  v2(39,1)=2;
  v2(39,2)=254;
  v2(39,3)=(-(T53*T47*(T149+T13*T222)));
  v2(40,1)=2;
  v2(40,2)=128;
  v2(40,3)=  v2(39,3);
  v2(41,1)=2;
  v2(41,2)=261;
  v2(41,3)=(-(T53*T47*(T13*T148+T13*T13*T221)));
  v2(42,1)=2;
  v2(42,2)=181;
  v2(42,3)=T42*params(3)*exp(y(10))+T39*params(3)*exp(y(10))+params(3)*(exp(y(10))-exp(params(4)))*T39+T39*params(3)*exp(y(10));
  v2(43,1)=2;
  v2(43,2)=271;
  v2(43,3)=(-(T133*T16*(T47+exp(y(15))*exp(y(15))*params(1)*params(3))));
  v2(44,1)=2;
  v2(44,2)=91;
  v2(44,3)=  v2(43,3);
  v2(45,1)=2;
  v2(45,2)=279;
  v2(45,3)=(-(T53*T16*(T47+exp(y(15))*exp(y(15))*params(1)*params(3))));
  v2(46,1)=2;
  v2(46,2)=243;
  v2(46,3)=  v2(45,3);
  v2(47,1)=2;
  v2(47,2)=273;
  v2(47,3)=(-(T53*T149*(T47+exp(y(15))*exp(y(15))*params(1)*params(3))));
  v2(48,1)=2;
  v2(48,2)=129;
  v2(48,3)=  v2(47,3);
  v2(49,1)=2;
  v2(49,2)=280;
  v2(49,3)=(-(T53*T13*T148*(T47+exp(y(15))*exp(y(15))*params(1)*params(3))));
  v2(50,1)=2;
  v2(50,2)=262;
  v2(50,3)=  v2(49,3);
  v2(51,1)=2;
  v2(51,2)=281;
  v2(51,3)=(-(T53*T16*(T47+exp(y(15))*exp(y(15))*params(1)*params(3)+exp(y(15))*exp(y(15))*params(1)*params(3)+exp(y(15))*exp(y(15))*params(1)*params(3))));
  v2(52,1)=3;
  v2(52,2)=21;
  v2(52,3)=(-(exp(x(it_, 1))*T66*(T169+exp(y(2))*exp(y(2))*getPowerDeriv(exp(y(2)),params(9),2))));
  v2(53,1)=3;
  v2(53,2)=141;
  v2(53,3)=exp(y(8));
  v2(54,1)=3;
  v2(54,2)=154;
  v2(54,3)=(-(exp(x(it_, 1))*T169*T174));
  v2(55,1)=3;
  v2(55,2)=28;
  v2(55,3)=  v2(54,3);
  v2(56,1)=3;
  v2(56,2)=161;
  v2(56,3)=(-(exp(x(it_, 1))*T69*(T174+exp(y(9))*exp(y(9))*getPowerDeriv(exp(y(9)),1-params(9),2))));
  v2(57,1)=3;
  v2(57,2)=306;
  v2(57,3)=(-(exp(x(it_, 1))*T66*T169));
  v2(58,1)=3;
  v2(58,2)=36;
  v2(58,3)=  v2(57,3);
  v2(59,1)=3;
  v2(59,2)=313;
  v2(59,3)=(-(exp(x(it_, 1))*T69*T174));
  v2(60,1)=3;
  v2(60,2)=169;
  v2(60,3)=  v2(59,3);
  v2(61,1)=3;
  v2(61,2)=321;
  v2(61,3)=(-(T66*T69*exp(x(it_, 1))));
  v2(62,1)=4;
  v2(62,2)=81;
  v2(62,3)=(-(T81*(T88*T136+T88*T315)));
  v2(63,1)=4;
  v2(63,2)=161;
  v2(63,3)=exp(y(9));
  v2(64,1)=4;
  v2(64,2)=176;
  v2(64,3)=(-(T88*T136*T184));
  v2(65,1)=4;
  v2(65,2)=86;
  v2(65,3)=  v2(64,3);
  v2(66,1)=4;
  v2(66,2)=181;
  v2(66,3)=(-(T90*T77*(T183+T78*T78*getPowerDeriv(T78,params(6),2))));
  v2(67,1)=4;
  v2(67,2)=195;
  v2(67,3)=(-(T81*(T206+T205*T315)));
  v2(68,1)=4;
  v2(68,2)=87;
  v2(68,3)=  v2(67,3);
  v2(69,1)=4;
  v2(69,2)=200;
  v2(69,3)=(-(T184*T206));
  v2(70,1)=4;
  v2(70,2)=182;
  v2(70,3)=  v2(69,3);
  v2(71,1)=4;
  v2(71,2)=201;
  v2(71,3)=(-(T81*(T205*T205*T314+T136*((-(exp(y(5))*T87))*T87*T87-(-(exp(y(5))*T87))*(T87*T87+T87*T87))/(T87*T87*T87*T87))));
  v2(72,1)=5;
  v2(72,2)=81;
  v2(72,3)=(exp(y(5))*exp(y(5))*(-(exp(y(7))*exp(y(5))))-(-(exp(y(7))*exp(y(5))))*(exp(y(5))*exp(y(5))+exp(y(5))*exp(y(5))))/(exp(y(5))*exp(y(5))*exp(y(5))*exp(y(5)));
  v2(73,1)=5;
  v2(73,2)=119;
  v2(73,3)=T142;
  v2(74,1)=5;
  v2(74,2)=83;
  v2(74,3)=  v2(73,3);
  v2(75,1)=5;
  v2(75,2)=121;
  v2(75,3)=T93;
  v2(76,1)=5;
  v2(76,2)=181;
  v2(76,3)=T97*(T78*2*(T78-1)+T78*2*T78);
  v2(77,1)=5;
  v2(77,2)=201;
  v2(77,3)=(-(((-exp(y(11)))*exp(y(11))*exp(y(11))-(-exp(y(11)))*(exp(y(11))*exp(y(11))+exp(y(11))*exp(y(11))))/(exp(y(11))*exp(y(11))*exp(y(11))*exp(y(11)))));
  v2(78,1)=8;
  v2(78,2)=1;
  v2(78,3)=(-exp(y(1)));
  v2(79,1)=8;
  v2(79,2)=101;
  v2(79,3)=exp(y(6));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),8,361);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,6859);
end
end
end
end

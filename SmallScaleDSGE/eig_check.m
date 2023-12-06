%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'eig_check';
M_.dynare_version = '4.5.7';
oo_.dynare_version = '4.5.7';
options_.dynare_version = '4.5.7';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = 'u_r';
M_.exo_names_tex = 'u\_r';
M_.exo_names_long = 'u_r';
M_.exo_names = char(M_.exo_names, 'u_g');
M_.exo_names_tex = char(M_.exo_names_tex, 'u\_g');
M_.exo_names_long = char(M_.exo_names_long, 'u_g');
M_.exo_names = char(M_.exo_names, 'u_z');
M_.exo_names_tex = char(M_.exo_names_tex, 'u\_z');
M_.exo_names_long = char(M_.exo_names_long, 'u_z');
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names_long = 'GDP';
M_.endo_names = char(M_.endo_names, 'x');
M_.endo_names_tex = char(M_.endo_names_tex, 'x');
M_.endo_names_long = char(M_.endo_names_long, 'x');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'Consumption');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'Nominal');
M_.endo_names = char(M_.endo_names, 'real');
M_.endo_names_tex = char(M_.endo_names_tex, 'real');
M_.endo_names_long = char(M_.endo_names_long, 'Real');
M_.endo_names = char(M_.endo_names, 'infl');
M_.endo_names_tex = char(M_.endo_names_tex, 'infl');
M_.endo_names_long = char(M_.endo_names_long, 'Inflation');
M_.endo_names = char(M_.endo_names, 'g');
M_.endo_names_tex = char(M_.endo_names_tex, 'g');
M_.endo_names_long = char(M_.endo_names_long, 'Fiscal');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'Deamnd');
M_.endo_partitions = struct();
M_.param_names = 'BETA';
M_.param_names_tex = 'BETA';
M_.param_names_long = 'BETA';
M_.param_names = char(M_.param_names, 'nu');
M_.param_names_tex = char(M_.param_names_tex, 'nu');
M_.param_names_long = char(M_.param_names_long, 'nu');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'piss');
M_.param_names_tex = char(M_.param_names_tex, 'piss');
M_.param_names_long = char(M_.param_names_long, 'piss');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, 'tau');
M_.param_names_long = char(M_.param_names_long, 'tau');
M_.param_names = char(M_.param_names, 'psi1');
M_.param_names_tex = char(M_.param_names_tex, 'psi1');
M_.param_names_long = char(M_.param_names_long, 'psi1');
M_.param_names = char(M_.param_names, 'psi2');
M_.param_names_tex = char(M_.param_names_tex, 'psi2');
M_.param_names_long = char(M_.param_names_long, 'psi2');
M_.param_names = char(M_.param_names, 'rhog');
M_.param_names_tex = char(M_.param_names_tex, 'rhog');
M_.param_names_long = char(M_.param_names_long, 'rhog');
M_.param_names = char(M_.param_names, 'rhoR');
M_.param_names_tex = char(M_.param_names_tex, 'rhoR');
M_.param_names_long = char(M_.param_names_long, 'rhoR');
M_.param_names = char(M_.param_names, 'gss');
M_.param_names_tex = char(M_.param_names_tex, 'gss');
M_.param_names_long = char(M_.param_names_long, 'gss');
M_.param_names = char(M_.param_names, 'yss');
M_.param_names_tex = char(M_.param_names_tex, 'yss');
M_.param_names_long = char(M_.param_names_long, 'yss');
M_.param_names = char(M_.param_names, 'gammaQ');
M_.param_names_tex = char(M_.param_names_tex, 'gammaQ');
M_.param_names_long = char(M_.param_names_long, 'gammaQ');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'zss');
M_.param_names_tex = char(M_.param_names_tex, 'zss');
M_.param_names_long = char(M_.param_names_long, 'zss');
M_.param_names = char(M_.param_names, 'zeta');
M_.param_names_tex = char(M_.param_names_tex, 'zeta');
M_.param_names_long = char(M_.param_names_long, 'zeta');
M_.param_names = char(M_.param_names, 'piA');
M_.param_names_tex = char(M_.param_names_tex, 'piA');
M_.param_names_long = char(M_.param_names_long, 'piA');
M_.param_names = char(M_.param_names, 'rhoz');
M_.param_names_tex = char(M_.param_names_tex, 'rhoz');
M_.param_names_long = char(M_.param_names_long, 'rhoz');
M_.param_names = char(M_.param_names, 'psi11');
M_.param_names_tex = char(M_.param_names_tex, 'psi11');
M_.param_names_long = char(M_.param_names_long, 'psi11');
M_.param_names = char(M_.param_names, 'psi12');
M_.param_names_tex = char(M_.param_names_tex, 'psi12');
M_.param_names_long = char(M_.param_names_long, 'psi12');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 8;
M_.param_nbr = 19;
M_.orig_endo_nbr = 8;
M_.aux_vars = [];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
erase_compiled_function('eig_check_static');
erase_compiled_function('eig_check_dynamic');
M_.orig_eq_nbr = 8;
M_.eq_nbr = 8;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 5 13;
 0 6 0;
 0 7 14;
 2 8 0;
 0 9 0;
 0 10 15;
 3 11 0;
 4 12 16;]';
M_.nstatic = 2;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 2;
M_.nsfwrd   = 4;
M_.nspred   = 4;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(19, 1);
M_.NNZDerivatives = [31; 79; -1];
load PARAMSTORE;
set_param_value('BETA',BETA);
set_param_value('nu',nu);
set_param_value('tau',tau);
set_param_value('rhoR',rhoR);
set_param_value('psi1',psi1);
set_param_value('psi2',psi2);
set_param_value('phi',phi);
set_param_value('rhog',rhog);
set_param_value('rhoz',rhoz);
set_param_value('zeta',zeta);
set_param_value('gammaQ',gammaQ);
set_param_value('piA',piA);
M_.params( 10 ) = log(1/(1-M_.params(15)));
gss = M_.params( 10 );
M_.params( 4 ) = log(1+M_.params(16)/400);
piss = M_.params( 4 );
M_.params( 11 ) = log((1-M_.params(2))^(1/M_.params(5))/(1-M_.params(15)));
yss = M_.params( 11 );
M_.params( 13 ) = 1+M_.params(12)/100;
gamma = M_.params( 13 );
M_.params( 14 ) = log(M_.params(13));
zss = M_.params( 14 );
r=log(gamma)-log(BETA);
M_.params( 18 ) = (1-M_.params(9))*M_.params(6);
psi11 = M_.params( 18 );
M_.params( 19 ) = (1-M_.params(9))*M_.params(7);
psi12 = M_.params( 19 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 3 ) = log((1-M_.params(2))^(1/M_.params(5)));
oo_.steady_state( 1 ) = M_.params(11);
oo_.steady_state( 5 ) = log(M_.params(13)/M_.params(1));
oo_.steady_state( 4 ) = log(M_.params(13)*exp(M_.params(4))/M_.params(1));
oo_.steady_state( 6 ) = log(1+M_.params(16)/400);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (sigma_r)^2;
M_.Sigma_e(2, 2) = (sigma_g)^2;
M_.Sigma_e(3, 3) = (sigma_z)^2;
options_.noprint=1;
oo_.dr.eigval = check(M_,options_,oo_);
save('eig_check_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('eig_check_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('eig_check_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('eig_check_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('eig_check_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('eig_check_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('eig_check_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end

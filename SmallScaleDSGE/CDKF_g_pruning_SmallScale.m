% This function evaluates the measurement equations 

function [y_bar] = CDKF_g_pruning_SmallScale(x_cu,A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt,state_var_positions)
    ne=size(S2,1);
    [~, ns] = size(B);
    x_f_cu=x_cu(1:ns,1);
    x_s_cu=x_cu(ns+1:end,1);
    y_bar=A+B*x_s_cu;

end
% % Setting the variable for reporting errors
% error_mes = 0;
% 
% % Determining the number of state variables to estimate
% xf_cu(:,1)     = x_cu(1:nx,1);
% xs_cu          = zeros(nx,1);
% xs_cu(1:nx1,1) = x_cu(nx+1:dimx,1);
% y              = zeros(dimy,1);
% AA             = 0.5*gss*sig^2;
% CC             = vech(xf_cu(1:nx,1)*xf_cu(1:nx,1)');
% for i=1:dimy
%     y(i,1) = gx(i,:)*(xf_cu+xs_cu) + GGxx(i,:)*CC + AA(i,1);     
% end

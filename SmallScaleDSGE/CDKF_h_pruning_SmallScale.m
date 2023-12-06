% This function evaluates the functions for the state equations 

function [x_cup] = CDKF_h_pruning_SmallScale(x_cu,v_cu,A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt,state_var_positions)
    ne=size(S2,1);
    [~, ns] = size(B);
    x_f_cu=x_cu(1:ns,1);
    x_s_cu=x_cu(ns+1:end,1);
    x_f=x_f_cu(state_var_positions,:);
    x_s=x_s_cu(state_var_positions,:);
    
    %the first order effects
    x_f_cup=Phi*x_f+R*v_cu;
    
    %the second order effects
    x_s_cup = 0.5*omega_ss_delta+Phi*x_s+R*v_cu+...
              0.5*Phi11*kron(x_f,x_f)+Phi12*kron(x_f,v_cu)+0.5*R00*kron(v_cu,v_cu);
    
    x_cup = [x_f_cup;x_s_cup];
    




end
% % Determining the number of state variables to estimate
% xf_cu(:,1)             = x_cu(1:nx,1);
% xs_cu(1:nx1,1)         = x_cu(nx+1:dimx,1);
% 
% % The first order effects
% xf_cup = hx*xf_cu;
% 
% % The second order effects
% xs_cup(1:nx1,1) = hx(1:nx1,1:nx1)*xs_cu(1:nx1,1) + HHxx(1:nx1,:)*vech(xf_cu*xf_cu') + 0.5*sig^2*hss(1:nx1,1);
%     
% % The state vector in the next period 
% x_cup(:,1) = [xf_cup(1:nx,1); xs_cup(1:nx1,1)];


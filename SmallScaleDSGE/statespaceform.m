function [omega_ss,F1,F2,omega_ss_delta,F11,F22,F12,F111,F112,F122,F222,F331,state_var_positions,control_var_postion] = statespaceform(optPruning,oo_,M_,f_11,order)

if order==1
    %First order solution is of the form:
    %\omega_t=\omega_ss+F1*\omega_{t-1}+F2*\epsilon_t
    %in dynare \omega_t=y;\epsilon_t=u;
    
    omega_ss=oo_.dr.ys;
    F1=oo_.dr.ghx;
    F2=oo_.dr.ghu;
    F11=0;
    F22=0;
    F12=0;
    F111=0;
    F112=0;
    F122=0;
    F222=0;
    F331=0;
   
    omega_ss_delta=0;
elseif order==2
    %Second order solution is of the form:
    %\omega_t=(\omega_ss+0.5*Delta^2)+F1*\omega_{t-1}+F2*\epsilon_t+0.5*F11*kron(x_{t-1},x_{t-1})+...
    %           0.5*F22*kron(\epsilon_t,\epsilon_t)+F12*kron(x_{t-1},\epsilon_t);
    %Delta^2=output.omega_ss_delta=adjustment to steady state due to variances of shocks
    omega_ss=oo_.dr.ys;
    omega_ss_delta=oo_.dr.ghs2;
    F1=oo_.dr.ghx;
    F2=oo_.dr.ghu;
    F11=oo_.dr.ghxx;
    F12=oo_.dr.ghxu;
    F22=oo_.dr.ghuu;
    F111=0;
    F112=0;
    F122=0;
    F222=0;
    F331=0;
   
elseif order==3
    [f_ss,fv,fvv,fss,fvvv,fssv,fsss,sigma,nv,nu,nx,nz] = RunDynarePruningpc(optPruning,oo_,M_,f_11);
    omega_ss=oo_.dr.ys;
    omega_ss_delta=oo_.dr.ghs2;
    F1=oo_.dr.ghx;
    F2=oo_.dr.ghu;
    F11=oo_.dr.ghxx;
    F12=oo_.dr.ghxu;
    F22=oo_.dr.ghuu;
    f111r=fvvv(:,1:nx,1:nx,1:nx);
    F111=reshape(f111r,M_.endo_nbr,nx^3);%div by 6
    f112r=fvvv(:,1:nx,1:nx,nx+1:nx+nu);
    F112=reshape(f112r,M_.endo_nbr,nx^2*nu);
    f122r=fvvv(:,1:nx,nx+1:nx+nu,nx+1:nx+nu);
    F122=reshape(f122r,M_.endo_nbr,nx*nu*nu);
    f222r=fvvv(:,nx+1:nx+nu,nx+1:nx+nu,nx+1:nx+nu);%div by 6
    F222=reshape(f222r,M_.endo_nbr,nu*nu*nu);
%     F331=fssv ;
    F331=6/3*(oo_.dr.g_1-[oo_.dr.ghx oo_.dr.ghu]);
end

% state_var_positions=oo_.dr.state_var;


for i=1:length(oo_.dr.order_var)
    for j=1:length(oo_.dr.state_var)
        if oo_.dr.state_var(1,j)==oo_.dr.order_var(i,1);
            state_var_positions(j,1)=(find((oo_.dr.order_var(:,1)==oo_.dr.state_var(1,j))));
        end
    end
end
control_var_postion=[];

% state_var_positions=[8,9,10];
% % allvars=1:length(M_.endo_names);
% % control_var_postion=[];
% % for i=1:length(allvars)
% %     if allvars(1,i)~=state_var_positions(1,:)
% %         control_var_postion=[control_var_postion,allvars(1,i)];
% %     end
% % end

end

% function [output,state_var_positions] = statespaceform(oo_,order)
%     output.omega_ss=oo_.dr.ys;
%     output.omega_ss_delta=oo_.dr.ghs2;
%     output.F1=oo_.dr.ghx;
%     output.F2=oo_.dr.ghu;
%     output.F11=oo_.dr.ghxx;
%     output.F12=oo_.dr.ghxu;
%     output.F22=oo_.dr.ghuu;

%     output.omega_ss=oo_.dr.ys;
%     output.F1=oo_.dr.ghx;
%     output.F2=oo_.dr.ghu;


function [xhat,Smat,LogL]=CDKF_dsge_DG_SmallScale(kalmfilex,kalmfiley,...
    A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt,state_var_positions)

%   This function performs the Central Difference Kalman filter; a state estimation for nonlinear 
%   systems that is based on second-order polynomial approximations of the 
%   nonlinear mappings. The approximations are derived by using a 
%   multidimensional extension of Stirling's interpolation formula. 
%   The model of the nonlinear system must be specified in the form:
%               x(k+1) = h[x(k)]+w(k+1)
%               y(k)   = g[x(k),z(k)]+v(k)
%   where 'x' is the state vector, y is the control vector
%   Finally, 'v' and 'w'  are (white) noise sources.
%
%
% Input
%  kalmfilex - File containing the state equations.
%  kalmfiley - File containing the output equations.
%  x0        - Initial state vector of dimension dimx*1
%  P0        - Initial covariance matrix (symmetric, nonnegative definite).
%  Rw        - Covariance matrices for w.
%  y         - Observed output. Dimension is [observations * outputs].
%  Rest      - Variables etc. related to the DSGE model
%
% Output
%   xhat    - State estimates. Dimension is [states, samples].
%             Notice: in xhat, the first row is for observation 0.
%   Smat    - Matrix with the square root of the covariance matrix for the
%             posterior state estimates. The dimension is [dimx, dimx, samples] 
%   LogL    - The contribution to SumLogL for each observation
%
%
%  The user must write the m-functions 'kalmfilex' and 'kalmfiley' 
%
%  Literature:
%     M. Norgaard, N.K. Poulsen, O. Ravn: "New Developments in State
%     Estimation for Nonlinear Systems", Automatica, (36:11), Nov. 2000,
%     pp. 1627-1638.
%
% These codes are modified by Martin Møller Andreasen (MMA) in spring 2008 based 
% on the codes for the dd2 filter by Norregaard et al. In this implementation, an 
% approximated expression for the quasi log-likelihood value is also calculated. 

% >>>>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<
h2        = 3;                        % Squared divided-difference step size
h         = sqrt(h2);                 % Divided difference step-size
scal1     = 0.5/h;                    % A convenient scaling factor
scal2     = sqrt((h2-1)/(4*h2*h2));   % Another scaling factor
dimy=size(yt,1);
T=size(yt,2);
%samples   = size(yt,1);                % The samples length
%dimy      = size(yt,2);                % # of observations at each time point
%dimx      = size(P0,1);               % # of states
ne=size(S2,1);
[~, ns] = size(B);
%[nx,ne]   = size(eta);                % # of state variables in the DSGE model and 
                                      % # of structural shocks in the DSGE model
%nx1       = nx-ne;                    % # of endogenous state variables in the DSGE model
%Rv        = zeros(dimx,dimx);         % The covariance matrix for the structural shocks
%for i=1:ne
%    Rv(nx1+i,nx1+i) = (sig*eta(nx1+i,i))^2;
%end
% x0_first  = zeros(ns,1);
% P0_first  = 10^(-6)*eye(ns);
% 
% x0_second  = zeros(ns,1);
% P0_second  = 10^(-6)*eye(ns);

x0        = zeros(2*ns,1);
P0        = 10^(-6)*eye(2*ns);
Rw        = S2;
dimv      = size(H,1);               % # of measurement noise sources
dimw      = size(Rw,1);               % # of structural shocks sources
[v,d]     = eig(P0);                  % Square root of initial state covariance
Sx        = triag(real(v*sqrt(d)));
[v,d]     = eig(Rw);                  % Square root of measurement error covariance
Sw        = real(v*sqrt(d));
[v,d]     = eig(H);                  % Square root of structural noise covariance
Sv        = real(v*sqrt(d));             

% Setting the dimensions
% SxxSxv    = zeros(2*ns,2*ns+2*ne);            % Allocate compund matrix consisting of Sxx and Sxv 
% SyxSyw    = zeros(dimy,2*ns+dimw);            % Allocate compund matrix consisting of Syx and Syw
xhat      = zeros(2*ns,T);                % Matrix for storing state estimates
Smat      = zeros(2*ns,2*ns,T);           % Matrix for storing cov. matrices
LogL      = zeros(T,1);                   % Vector for storing the contributions to the log likelihood fct

% Linear process noise model in state equation: x(k+1) = f(x(k))+w(k+1)


% Linear observation noise model: y(k) = g(x(k),z(k)) + v(k)
%SyxSyw(:,ns+1:ns+dimv) = Sv;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>> FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<
for t=1:T
    if t == 1
       fxbar = feval(kalmfilex,x0,zeros(ne,1),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
       xbar = ((h2-2*ns-ne)/h2)*fxbar;
       for kx=1:2*ns
           sxp              = feval(kalmfilex,x0+h*Sx(:,kx),zeros(ne,1),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
           sxm              = feval(kalmfilex,x0-h*Sx(:,kx),zeros(ne,1),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
           Sxx_first(:,kx)  = scal1*(sxp-sxm);
           Sxx_second(:,kx) = scal2*(sxp+sxm-2*fxbar); 
           xbar             = xbar + (sxp+sxm)/(2*h2);
       end
       
       for kx=1:ne
           swp              = feval(kalmfilex,x0,zeros(ne,1)+h*Sw(:,kx),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
           swm              = feval(kalmfilex,x0,zeros(ne,1)-h*Sw(:,kx),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
           Sxw_first(:,kx)  = scal1*(swp-swm);
           Sxw_second(:,kx) = scal2*(swp+swm-2*fxbar);
           xbar             = xbar + (swp+swm)/(2*h2);
       end
       
       
    else
        fxbar  = feval(kalmfilex,xhat(:,t-1),zeros(ne,1),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);    
        xbar = ((h2-2*ns-ne)/h2)*fxbar;
        for kx=1:2*ns
            sxp              = feval(kalmfilex,xhat(:,t-1)+h*Sx(:,kx),zeros(ne,1),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
            sxm              = feval(kalmfilex,xhat(:,t-1)-h*Sx(:,kx),zeros(ne,1),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
            Sxx_first(:,kx)  = scal1*(sxp-sxm);
            Sxx_second(:,kx) = scal2*(sxp+sxm-2*fxbar); 
            xbar             = xbar + (sxp+sxm)/(2*h2);
        end
        
        for kx=1:ne
           swp              = feval(kalmfilex,xhat(:,t-1),zeros(ne,1)+h*Sw(:,kx),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
           swm              = feval(kalmfilex,xhat(:,t-1),zeros(ne,1)-h*Sw(:,kx),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
           Sxw_first(:,kx)  = scal1*(swp-swm);
           Sxw_second(:,kx) = scal2*(swp+swm-2*fxbar);
           xbar             = xbar + (swp+swm)/(2*h2);
       end
        
    end     
%    % Cholesky factor of a'priori estimation error covariance
     Sxbar = triag([Sxx_first,Sxw_first,Sxx_second,Sxw_second]);
%     
%     % --- Measurement update (a posteriori update) ---
     [y0] = feval(kalmfiley,xbar,A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
     ybar = ((h2-2*ns)/h2)*y0;
%    kx2 = dimx+dimw;
     for kx=1:2*ns
        [syp]            = feval(kalmfiley,xbar+h*Sxbar(:,kx),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
        [sym]            = feval(kalmfiley,xbar-h*Sxbar(:,kx),A,B,H,Phi,R,S2,Phi11,Phi12,R00,omega_ss_delta,yt(:,t),state_var_positions);
        Syx_first(:,kx)  = scal1*(syp-sym);
        Syx_second(:,kx) = scal2*(syp+sym-2*y0);
        ybar             = ybar+(syp+sym)/(2*h2);    
     end
% 
%    % Cholesky factor of a'posteriori output estimation error covariance
     Sy      = triag([Syx_first,Sv,Syx_second]);
     Pyy     = Sy*Sy';
     inv_Pyy = Pyy\eye(dimy);
     Kgain   = (Sxbar*Syx_first')*inv_Pyy;
%          
%    % Contribution to the quasi log-likelihood function
     LogL_indv(t,1) = -0.5*dimy*log(2*pi)-0.5*log(det(Pyy))...
           -0.5*(yt(:,t)'-ybar(:,1)')*inv_Pyy*(yt(:,t)-ybar(:,1));
%                
%    % Update of the state estimate - i.e. posterior state estimate
     xhat(:,t) = xbar + Kgain*(yt(:,t)-ybar(:,1)); 
         
%    % Cholesky factor of a'posteriori estimation error covariance
     Sx    = triag([Sxbar-Kgain*Syx_first,Kgain*Sv,Kgain*Syx_second]); 
     %Sx   = triag([Sxbar-Kgain*SyxSyw(:,1:dimx) Kgain*SyxSyw(:,dimx+1:end)]);
% 
%    % Store the covariances for the state estimate
     Smat(:,:,t) = Sx;
end

LogL=sum(LogL_indv);
end
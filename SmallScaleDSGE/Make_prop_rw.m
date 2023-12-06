load('PMMH_SmallScale_v4_025.mat');
burn=10000;
theta=[log(Post.tau(burn:end,1)),log(Post.psi1(burn:end,1)),log(Post.psi2(burn:end,1)),logit_inverse(Post.rhoR(burn:end,1)),...
    logit_inverse(Post.rhog(burn:end,1)),logit_inverse(Post.rhoz(burn:end,1)),log(Post.sigma_r(burn:end,1)),log(Post.sigma_g(burn:end,1)),log(Post.sigma_z(burn:end,1)),...
    log(Post.rA(burn:end,1)),log(Post.piA(burn:end,1)),Post.gammaQ(burn:end,1)];
thetasave=theta;
covmat=cov(thetasave);
covmat=jitChol(covmat);
save('covmat_smallscale_2ndorder_025.mat','covmat');



load('PMMH_SmallScale_v4_nonfixed.mat');
burn=19000;
theta=[log(Post.tau(burn:end,1)),log(Post.psi1(burn:end,1)),log(Post.psi2(burn:end,1)),logit_inverse(Post.rhoR(burn:end,1)),...
    logit_inverse(Post.rhog(burn:end,1)),logit_inverse(Post.rhoz(burn:end,1)),log(Post.sigma_r(burn:end,1)),log(Post.sigma_g(burn:end,1)),log(Post.sigma_z(burn:end,1)),...
    log(Post.rA(burn:end,1)),log(Post.piA(burn:end,1)),Post.gammaQ(burn:end,1),...
    log(Post.sigma_measure1(burn:end,1)),log(Post.sigma_measure2(burn:end,1)),log(Post.sigma_measure3(burn:end,1))];
thetasave=theta;
covmat=cov(thetasave);
covmat=jitChol(covmat);
save('covmat_smallscale_2ndorder_nonfixed.mat','covmat');



%theta=[log(param_tau),log(param_psi1),log(param_psi2),logit_inverse(param_rhoR),...
%       logit_inverse(param_rhog),logit_inverse(param_rhoz),log(param_sigma_r),log(param_sigma_g),log(param_sigma_z),...
%       log(param_rA),log(param_piA),gammaQ,log(param_sigma_measure1),log(param_sigma_measure2),log(param_sigma_measure3)];
    
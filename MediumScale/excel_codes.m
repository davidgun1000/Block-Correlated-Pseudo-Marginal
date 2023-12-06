
% load('PMMH_2ndorder_05_T200_N100_optscale.mat');
% xlswrite('MPM_T200_N100_optscale.xlsx',[Post.RHO(5001:end,1),Post.RHOE(5001:end,1),Post.sigma_a(5001:end,1),Post.sigma_e(5001:end,1)]);
% load('PMMH_2ndorder_05_T200_N250_optscale.mat');
% xlswrite('MPM_T200_N250_optscale.xlsx',[Post.RHO(2001:end,1),Post.RHOE(2001:end,1),Post.sigma_a(2001:end,1),Post.sigma_e(2001:end,1)]);
% load('corrPMMH_2ndorder_05_T200_N100_optscale.mat');
% xlswrite('corrPMMH_N100_optscale.xlsx',[Post.RHO(5001:end,1),Post.RHOE(5001:end,1),Post.sigma_a(5001:end,1),Post.sigma_e(5001:end,1)]);
% load('corrPMMH_2ndorder_05_T200_N1000_optscale.mat');
% xlswrite('corrPMMH_N1000_optscale.xlsx',[Post.RHO(5001:end,1),Post.RHOE(5001:end,1),Post.sigma_a(5001:end,1),Post.sigma_e(5001:end,1)]);
% 
% xlswrite('PMMH_SmallScale_MeaError_prop1_N100_2_NIASRA.xlsx',[Post.piA,Post.tau,Post.psi1,Post.psi2,Post.gammaQ,Post.rA,Post.rhoR,Post.rhog,Post.rhoz,Post.sigma_r,Post.sigma_g,Post.sigma_z,Post.sigma_measure1,Post.sigma_measure2,...
%     Post.sigma_measure3])
% 
% xlswrite('PMMH_SmallScale_MeaError_prop1_N100_2.xlsx',[Post.piA(1:80000,1),Post.tau(1:80000,1),Post.psi1(1:80000,1),Post.psi2(1:80000,1),Post.gammaQ(1:80000,1),Post.rA(1:80000,1),Post.rhoR(1:80000,1),Post.rhog(1:80000,1),Post.rhoz(1:80000,1),Post.sigma_r(1:80000,1),Post.sigma_g(1:80000,1),Post.sigma_z(1:80000,1),Post.sigma_measure1(1:80000,1),Post.sigma_measure2(1:80000,1),...
%     Post.sigma_measure3(1:80000,1)])
% 
% xlswrite('PMMH_SmallScale_MeaError_prop1_N250_NIASRA.xlsx',[Post.piA(7501:end,1),Post.tau(7501:end,1),Post.psi1(7501:end,1),Post.psi2(7501:end,1),Post.gammaQ(7501:end,1),Post.rA(7501:end,1),Post.rhoR(7501:end,1),Post.rhog(7501:end,1),Post.rhoz(7501:end,1),Post.sigma_r(7501:end,1),Post.sigma_g(7501:end,1),Post.sigma_z(7501:end,1),Post.sigma_measure1(7501:end,1),Post.sigma_measure2(7501:end,1),...
%     Post.sigma_measure3(7501:end,1)])
% 
% save('PMMH_SmallScale_MeaError_prop1_N250_NIASRA_Pratiti.mat',Post.piA(7501:end,1),Post.tau(7501:end,1),Post.psi1(7501:end,1),Post.psi2(7501:end,1),Post.gammaQ(7501:end,1),Post.rA(7501:end,1),Post.rhoR(7501:end,1),Post.rhog(7501:end,1),Post.rhoz(7501:end,1),Post.sigma_r(7501:end,1),Post.sigma_g(7501:end,1),Post.sigma_z(7501:end,1),Post.sigma_measure1(7501:end,1),Post.sigma_measure2(7501:end,1),...
%     Post.sigma_measure3(7501:end,1))

xlswrite('PMMH_SmallScale_v4_05_USED.xlsx',[Post.piA(:,1),Post.tau(:,1),Post.psi1(:,1),Post.psi2(:,1),Post.gammaQ(:,1),Post.rA(:,1),Post.rhoR(:,1),Post.rhog(:,1),Post.rhoz(:,1),Post.sigma_r(:,1),Post.sigma_g(:,1),Post.sigma_z(:,1),Post.sigma_measure1(:,1),Post.sigma_measure2(:,1),...
     Post.sigma_measure3(:,1)])

xlswrite('PMMH_SmallScale_v4_025_USED.xlsx',[Post.piA(:,1),Post.tau(:,1),Post.psi1(:,1),Post.psi2(:,1),Post.gammaQ(:,1),Post.rA(:,1),Post.rhoR(:,1),Post.rhog(:,1),Post.rhoz(:,1),Post.sigma_r(:,1),Post.sigma_g(:,1),Post.sigma_z(:,1),Post.sigma_measure1(:,1),Post.sigma_measure2(:,1),...
     Post.sigma_measure3(:,1)])

xlswrite('PMMH_SmallScale_v4_nonfixed_USED.xlsx',[Post.piA(:,1),Post.tau(:,1),Post.psi1(:,1),Post.psi2(:,1),Post.gammaQ(:,1),Post.rA(:,1),Post.rhoR(:,1),Post.rhog(:,1),Post.rhoz(:,1),Post.sigma_r(:,1),Post.sigma_g(:,1),Post.sigma_z(:,1),Post.sigma_measure1(:,1),Post.sigma_measure2(:,1),...
     Post.sigma_measure3(:,1)])

 




% xlswrite('sim_S19_T1000_mu.xlsx',theta_mu_store)
% xlswrite('sim_S19_T1000_sig2.xlsx',[theta_sig2_store1,theta_sig2_store2,theta_sig2_store3,theta_sig2_store4,theta_sig2_store5,theta_sig2_store6,theta_sig2_store7]);
% xlswrite('sim_S19_T1000_b1.xlsx',[theta_latent_b1_store]);
% xlswrite('sim_S19_T1000_b2.xlsx',[theta_latent_b2_store]);
% xlswrite('sim_S19_T1000_b3.xlsx',[theta_latent_b3_store]);
% xlswrite('sim_S19_T1000_A.xlsx',[theta_latent_A_store]);
% xlswrite('sim_S19_T1000_v1.xlsx',[theta_latent_v1_store]);
% xlswrite('sim_S19_T1000_v2.xlsx',[theta_latent_v2_store]);
% xlswrite('sim_S19_T1000_tau.xlsx',[theta_latent_tau_store]);


% xlswrite('PG_real_leverage_N100_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% xlswrite('PG_real_leverage_N250_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% xlswrite('PG_real_leverage_N500_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% xlswrite('PG_real_leverage_N1000_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% 
% xlswrite('mix_real_leverage_N100_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% xlswrite('mix_real_leverage_N250_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% xlswrite('mix_real_leverage_N500_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% xlswrite('mix_real_leverage_N1000_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% 
% xlswrite('corrmix_real_leverage_N100_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% 
% xlswrite('PG_Fearnhead_real_leverage_N100_T1000.xlsx',[Post.B1(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%     Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);
% xlswrite('PG_Fearnhead_real_leverage_N250_T1000.xlsx',[Post.B1(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%     Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);
% xlswrite('PG_Fearnhead_real_leverage_N500_T1000.xlsx',[Post.B1(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%     Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);
% xlswrite('PG_Fearnhead_real_leverage_N1000_T1000.xlsx',[Post.B1(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%     Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);
% 
% xlswrite('PG_Fearnhead_real_leverage_N500_T3000.xlsx',[Post.B1(1001:end,:),Post.B2(1001:end,:),Post.B3(1001:end,:),Post.B4(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%      Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);
% xlswrite('PG_Fearnhead_real_leverage_N1000_T3000.xlsx',[Post.B1(1001:end,:),Post.B2(1001:end,:),Post.B3(1001:end,:),Post.B4(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%      Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);
% xlswrite('PG_Fearnhead_real_leverage_N2000_T3000.xlsx',[Post.B1(1001:end,:),Post.B2(1001:end,:),Post.B3(1001:end,:),Post.B4(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%      Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);
% xlswrite('PG_Fearnhead_real_leverage_N5000_T3000.xlsx',[Post.B1(1001:end,:),Post.B2(1001:end,:),Post.B3(1001:end,:),Post.B4(1001:end,:),Post.mu(1001:end,:),Post.tau(1001:end,:),Post.phi(1001:end,:),Post.rho(1001:end,:),Post.ctraj1(1001:end,:),Post.ctraj2(1001:end,:),Post.ctraj10(1001:end,:),...
%      Post.ctraj11(1001:end,:),Post.ctraj12(1001:end,:),Post.ctrajfactor(1001:end,:)]);



% xlswrite('corr_real_leverage_N20.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_leverage_N50.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_leverage_N100.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_leverage_N200.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% 
% xlswrite('PG_real_leverage_N20.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_leverage_N50.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_leverage_N100.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_leverage_N200.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% 
% xlswrite('corr_real_leverage_N20_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('corr_real_leverage_N50_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('corr_real_leverage_N100_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('corr_real_leverage_N200_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% 
% xlswrite('PG_real_leverage_N20_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('PG_real_leverage_N50_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('PG_real_leverage_N100_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('PG_real_leverage_N200_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% 


% xlswrite('mix_real_N100.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('mix_real_N250.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('mix_real_N500.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('mix_real_N1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_EIS_real_N100.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_N100_leverage.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_N100_leverage.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);


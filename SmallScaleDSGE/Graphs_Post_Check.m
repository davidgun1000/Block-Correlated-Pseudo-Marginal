
load PMMH_SmallScale_MeaError_00020003_0502

draws=size(Post.sigma_r,1)
burn=0.20*draws

figure
subplot(3,1,1)
plot(Post.sigma_r(burn+1:end))
axis('tight')
subplot(3,1,2)
plot(Post.sigma_g(burn+1:end))
axis('tight')
subplot(3,1,3)
plot(Post.sigma_z(burn+1:end))
axis('tight')

mean(Post.sigma_r(burn+1:end))
mean(Post.sigma_g(burn+1:end))
mean(Post.sigma_z(burn+1:end))



figure
subplot(3,1,1)
plot(Post.rhoR(burn+1:end))
axis('tight')
subplot(3,1,2)
plot(Post.rhog(burn+1:end))
axis('tight')
subplot(3,1,3)
plot(Post.rhoz(burn+1:end))
axis('tight')


mean(Post.rhoR(burn+1:end))
mean(Post.rhog(burn+1:end))
mean(Post.rhoz(burn+1:end))

figure
subplot(3,1,1)
plot(Post.tau(burn+1:end))
axis('tight')
subplot(3,1,2)
plot(Post.psi1(burn+1:end))
axis('tight')
subplot(3,1,3)
plot(Post.psi2(burn+1:end))
axis('tight')


mean(Post.tau(burn+1:end))
mean(Post.psi1(burn+1:end))
mean(Post.psi2(burn+1:end))

figure
ksdensity(Post.rhog(burn+1:end))

mean(Post.sigma_measure1(burn+1:end))
mean(Post.sigma_measure2(burn+1:end))
mean(Post.sigma_measure3(burn+1:end))

%chheck the priors and posterior





mean(Post.tau(20001:end))

mean(Post.psi1(20001:end))

mean(Post.psi2(20001:end))






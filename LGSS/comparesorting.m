%comparing sorting method

N=2000;
d=30;

particles = randn(d,N);
weight = randn(1,N);
weight = weight./sum(weight);
u_res = rand(1,N);

for t=1:10
    t
tic
[indx]=rs_multinomial_corr2_PMMH_2ndorder_simplified(particles,weight,u_res);
time_DG(t,1)=toc;
end


for t=1:10
    t
tic
[indx]=rs_multinomial_corr2_PMMH_2ndorder(particles,weight,u_res);
time_C(t,1)=toc;
end

for t=1:10
    t
tic
[indx]=rs_multinomial_corr2_PMMH_2ndorder_simplified_Hilbert(particles,weight,u_res);
time_H(t,1)=toc;
end
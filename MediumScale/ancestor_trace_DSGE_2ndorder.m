function [sirhat_eps,sirhat_state_first,sirhat_state_second]=ancestor_trace_DSGE_2ndorder(particles_eps,particles1,particles2,w,indx)
%keyboard
% BACKTRACING (Kitagawa 1996, Andrieu et al, 2010)
T=size(indx,1);
N=size(indx,2);
outndx=NaN(1,T);
outndx(T)=randsample(N,1,true,w(T,:));
sirhat_state_first(:,T)=particles1(:,outndx(T),T);
sirhat_state_second(:,T)=particles2(:,outndx(T),T);
sirhat_eps(:,T)=particles_eps(:,outndx(T),T);
for t=T-1:-1:1,
    outndx(t)=indx(t+1,outndx(t+1));
    sirhat_state_first(:,t)=particles1(:,outndx(t),t);
    sirhat_state_second(:,t)=particles2(:,outndx(t),t);
    sirhat_eps(:,t)=particles_eps(:,outndx(t),t);
end
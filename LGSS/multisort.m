function [newp,neww,neworig]=multisort(vals,particles,w,indexorig)
nds=hilbertsort(vals); % sorting indices
newp=particles(:,nds); % sort particles
neww=w(:,nds); % sort weights
neworig=indexorig(:,nds);


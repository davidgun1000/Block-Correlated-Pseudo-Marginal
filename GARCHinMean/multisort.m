function [newp,neww,new_origindex]=multisort(particles,w,orig_index)
nds=hilbertsort(particles); % sorting indices
newp=particles(:,nds); % sort particles
neww=w(:,nds); % sort weights
new_origindex=orig_index(:,nds);
end
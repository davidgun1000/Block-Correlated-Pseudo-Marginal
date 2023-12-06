function [indx]=rs_multinomial_corr2_PMMH_2ndorder_simplified(particles,w,u)
N=length(w);
orig_index=(1:1:N);
[particles_sort,weight_sort,orig_ind_sort]=multidimsort2_simplified(particles,w,orig_index);
Q=cumsum(weight_sort); % cumulative sum
Q(1,N)=1;
u=sort(u);

%[~,~,indx_sort]=histcounts(u',[0,Q]');
%[~,indx_sort]=histc(u',cumsum([0;weight_sort']));


j=1;
for i=1:N
    while (Q(j)<u(i))
        j=j+1; % climb the ladder
    end
    indx_sort(i)=j; % assign index
end
indx=orig_ind_sort(indx_sort');
indx=indx';


end

%particles_sort=particles_sort';

% i=1;
% while (i<=N),
%     %sampl = rand(1,1);  % (0,1]
%     j=1;
%     while (Q(j)<u(i)),
%         j=j+1;
%     end;
%     indx_sort(i)=j;
%     i=i+1;
% end


%for i=1:N
%    tol=eps(particles_sort(indx_sort(i)));
%    id1=ismembertol(particles',particles_sort(indx_sort(i))*ones(N,1),tol);
%    indx(1,i)=find(id1>0);
%end

%for i=1:N
%    particles(1,i)=particles_sort(indx(i),1);
%end



% if conditional
%     indx(1)=1; % static conditional index
% end

% j=1;
% for i=1:N
%   while (Q(j)<u(i,1))
%     j=j+1;
%   end
%   indx(i)=j;
% end
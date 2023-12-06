function y=logitcdf_upperlower(x,lower,upper)
% input x in (-inf,inf)
%global PhiMax
y=1./(1+exp(-x));
y=lower+(upper-lower).*y;


end
%y=1./(1+exp(-x)); % transform in (0,1)
%y=PhiMax*y; % transform in (0,phimax)


function y=logitinverse_upperlower(x,lower,upper)
% input x in (0,phimax)

x=(x-lower)./(upper-lower);
y=log(x./(1-x));


end
%x=x/PhiMax; % transform in (0,1)
%y=log(x./(1-x)); % transform to (-inf,inf)

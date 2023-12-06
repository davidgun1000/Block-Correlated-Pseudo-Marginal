function [logpdf]=log_truncatednormal(x,mu,std_ev,lo,up)

    if x>lo & x<up
    x_std=(x-mu)./std_ev;
    lo_std=(lo-mu)./std_ev;
    up_std=(up-mu)./std_ev;
    
    num=log(normpdf(x_std));
    den=log(std_ev.*(normcdf(up_std)-normcdf(lo_std)));
    logpdf=num-den;
    else
    logpdf=-Inf;
    
    end




end
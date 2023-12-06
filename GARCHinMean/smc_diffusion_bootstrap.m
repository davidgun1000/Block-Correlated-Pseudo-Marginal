function [particles_eps,w,indx,sir_llh]=smc_diffusion_bootstrap(y,T,M,N,delta,a,tau,mu,A,u_particles,u_res,corr_index,disturbance_index)

   t=1;
   h=delta/M;
   ne=size(y,1);
   dim_y=size(y,1);
   particles_eps=zeros(ne,N,(M)*(T-1)+1);
   w=zeros((M)*(T-1)+1,N);
   indx=zeros((M)*(T-1)+1,N);
   particles_eps(:,1:N,t)=u_particles(:,1:N,t);
   particles_state(:,1:N,t)=particles_eps(:,1:N,t);
   r=1;
   error = (y(:,r) - A*particles_state(:,1:N,t));
   for i=1:dim_y
       logw_indv(i,:)=-0.5*log(2*pi)-0.5.*log(exp(particles_state(i,:,t)))-0.5.*(1./exp(particles_state(i,:,t))).*(error(i,:).^2);
   end
   logw=sum(logw_indv);
   w(t,:)=exp(logw-max(logw));
   w(t,:)=real(w(t,:));
   sir_llh(t,1)=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:));
   ESS=1/sum(w(t,:).^2);
   
   resample = 0;
   for t=2:((M)*(T-1)+1)
        
        if mod(t,M)==2 
           %[indx(t,:)]=[1:1:N]'; 
           resample = resample+1;
           if corr_index==1
              if disturbance_index==1    
                 particles_com=[particles_eps(:,1:N,t-1)];
              else
                 particles_com=[particles_state(:,1:N,t-1)];
              end 
              [indx_temp(t,:)]=rs_multinomial_corr2_PMMH_2ndorder_simplified(particles_com,w(t-1,:),u_res(t,:));
    
           else
           [indx_temp(t,:)]=rs_multinomial(w(t-1,:));
           end
           indx(t,:)=[indx_temp(t,:)];
           w(t-1,:) = 1/N*(ones(1,N));
        else
            [indx(t,:)]=[1:1:N]';
        end
        particles_res = particles_state(:,indx(t,1:N),t-1);
        z=particles_res;
%        %GARCH
        particles_eps(:,1:N,t)=u_particles(:,1:N,t);
        z=z+h.*(a.*(mu-exp(z)).*exp(-z)-(tau/2))+sqrt(tau)*sqrt(h)*particles_eps(:,1:N,t);
        z=real(z);
        particles_state(:,1:N,t)=z;
        
        if mod(t,M)==1
           r=r+1; 
           error = (y(:,r) - A*particles_state(:,1:N,t));
           for i=1:dim_y
               logw_indv(i,:)=-0.5*log(2*pi)-0.5.*log(exp(particles_state(i,:,t)))-0.5.*(1./exp(particles_state(i,:,t))).*(error(i,:).^2);
           end
           logw=sum(logw_indv);
           
           w(t,:)=w(t-1,:).*exp(logw-max(logw));
           %sir_llh=sir_llh+log(sum(w(t,:)))+max(logw);
           sir_llh(r,1)=log(sum(w(t,:)))+max(logw);
           
           w(t,:)=w(t,:)/sum(w(t,:));
           ESS=1/sum(w(t,:).^2);
        else
           w(t,:)=w(t-1,:);
           ESS=1/sum(w(t,:).^2);
        end
   end

end
%           logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,r)-B(1,:)*ft(:,r)).^2);


%    for t=2:T
%        if corr_index==1
%           if disturbance_index==1    
%              particles_com=[particles_eps(:,1:N,t-1)];
%           else
%              particles_com=[particles_state(:,1:N,t-1)];
%           end 
%           [indx_temp(t,:)]=rs_multinomial_corr2_PMMH_2ndorder_simplified(particles_com,w(t-1,:),u_res(t,:));
%        else
%        [indx_temp(t,:)]=rs_multinomial(w(t-1,:));
%        end
%        indx(t,:)=[indx_temp(t,:)]; 
%        particles_res=particles(:,indx(t,1:N),t-1);
%        z=particles_res;
%        for i=1:M
%            %OU model
%            %z=z+h*(a*(mu-z))+sqrt(tau)*sqrt(h)*randn(1,N);
%            
%            %GARCH diffusion model
%            z=z+h.*(a.*(mu-exp(z)).*exp(-z)-(tau/2))+sqrt(tau)*sqrt(h)*randn(1,N);
%            z=real(z);
%        end
%        particles(:,1:N,t)=z;
%        logw=-0.5*log(2*pi)-0.5.*log((exp(particles(1,:,t))))-0.5.*(1./(exp(particles(1,:,t)))).*((y(1,t)-B(1,:)*ft(:,t)).^2);
%        w(t,:)=exp(logw-max(logw));
%        sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
%        w(t,:)=w(t,:)/sum(w(t,:)); 
%    end
  







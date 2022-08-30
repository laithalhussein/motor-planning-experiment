function y = s2_dropped_nonlin_PO_model(b,x)
%X is passed in as:
%X = [m1,s1,s2, s1*0+st_std_gm, s2*0+dt_std_gm];
%SD ratio based on population means (all elements are the same value)
%mu2i = 0.5 * ( K1*(s2i) + (1-K1)*mean(s2) ) * ( K2*mu1i + (1-K2)*mean(mu1) )*( K3*(1/s1i) + (1-K3)*mean(1/s1) ) + K4

%Notes: 
%1) Make 0.5 factor a parameter
%2)Dont worry about *I* immediately, but if I is used, we should have a version where I is
%   based on a *scaled* SM and non-scaled
%3) Try when phi is subtracted from SM and without

mu1i = x(:,1);
mean_mu1 = mean(x(:,1));

s1i_inv = 1./x(:,2);
mean_s1_inv = mean(1./x(:,4));
%mean_s1_inv = mean(1./x(:,2));

s2i = x(:,3);
mean_s2 = mean(x(:,4));
%mean_s2 = mean(x(:,2));

obs_locy = 0.5 * 20 * sind(60);
phi = atand(2/(obs_locy - sind(60)));
I = abs(mu1i)>=phi;

mu1i = mu1i-phi;
mean_mu1 = mean_mu1-phi;

%include effect of I?
I = ones(size(I));

y = I.*( (0*s2i + mean_s2) .* (b(1)*mu1i + (1-b(1))*mean_mu1) .* (b(2)*s1i_inv + (1-b(2))*mean_s1_inv) )* b(3) + b(4);
end
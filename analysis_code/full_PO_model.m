function y = full_PO_model(b,x)
%first column is 1-tgt trial IMD, 2nd column is SD ratio, 3rd column is the
%SD ratio based on population means (all elements are the same value)
%mu2i = 0.5 * ( K1*(s2i) + (1-K1)*mean(s2) ) * ( K2*mu1i + (1-K2)*mean(mu1) )*( K3*(1/s1i) + (1-K3)*mean(1/s1) ) + K4

%Notes:Make 0.5 factor a parameter
%Dont worry about *I* immediately, but we should have a version where I is
%based on *scaled* SM

obs_locy = 0.5 * 20 * sind(60);
phi = atand(2/(obs_locy - sind(60)));
I = abs(x(:,1))>=phi;

x(:,1) = x(:,1)-phi;

I = ones(size(I));
y = 0.5*I.*( (mean(x(:,1)) + (1-b(1)).*(mean(x(:,1))-x(:,1))) .* (mean(x(:,3)) + (1-b(2)).*(mean(x(:,3))-x(:,2))) ) + b(3);
end
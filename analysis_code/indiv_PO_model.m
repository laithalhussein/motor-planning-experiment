function y = indiv_PO_model(b,x)
%first column is 1-tgt trial IMD, 2nd column is SD ratio, 3rd column is the
%SD ratio based on population means (all elements are the same value)
%m2_i = 0.5*( mean(m1) + alpha*(m1_i-mean(m1) ) * (mean(s) + beta*(s_i-mean(s) ) + offset
obs_locy = 0.5 * 20 * sind(60);
phi = atand(2/(obs_locy - sind(60)));
I = abs(x(:,1))>=phi;

x(:,1) = x(:,1)-phi;

I = ones(size(I));
y = 0.5*I.*x(:,1) .* x(:,2) + b(1);
end
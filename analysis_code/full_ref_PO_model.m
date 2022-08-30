function y = full_ref_PO_model(b,x)
%first column is 1-tgt trial IMD, 2nd column is SD ratio
%m2_i = K*0.5*( mean(m1) + alpha*(m1_i-mean(m1) ) * (mean(s) + beta*(s_i-mean(s) ) + offset
obs_locy = 0.5 * 20 * sind(60);
phi = atand(2/(obs_locy - sind(60)));
I = abs(x(:,1))>=phi;
I = ones(size(I));
y = b(1)*0.5*I.*( (mean(x(:,1)) + (1-b(2)).*(mean(x(:,1))-x(:,1))) .* (mean(x(:,2)) + (1-b(3)).*(mean(x(:,2))-x(:,2)))  + b(4) );
end
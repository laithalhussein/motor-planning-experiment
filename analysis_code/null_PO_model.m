function y = null_PO_model(b,x)
obs_locy = 0.5 * 20 * sind(60);
phi = atand(2/(obs_locy - sind(60)));
I = abs(x(:,1))>=phi;
I = ones(size(I));
x(:,1) = x(:,1)-phi;
y = 0.5*I.*( (mean(x(:,1)) + 0*(mean(x(:,1))-x(:,1))) .* (mean(x(:,3)) + 0.*(mean(x(:,3))-x(:,2))) + b);

end
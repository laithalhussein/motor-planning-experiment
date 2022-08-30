function y = single_dof_po_model3(b,x)
%first column is 1-tgt trial IMD, 2nd column is SD ratio
obs_locy = 0.5 * 20 * sind(60);
phi = atand(2/(obs_locy - sind(60)));
x(:,1) = -x(:,1);
I = abs(x(:,1))>=phi;
%I = ones(size(I));

x(:,1)= x(:,1)-phi;
y = b*I.*( x(:,1));
end
function y = full_mv1_PO_model(b,x)
%first column is 1-tgt trial IMD, 2nd column is SD on 1 tgt trials, 3rd
%column is SD on 2 tgt trials, 4th column is population-averaged estimate
%of SD on 1-target trials, and 5th column is population-averaged estimate
%of SD on 2-target trials
%%m2i = 0.5*( m1*s2 + alpha*( m1*s2 - m1i*s2i ) * ( 1/s1 + beta*(1/s1 - 1/s1i) ) + offset
obs_locy = 0.5 * 20 * sind(60);
phi = atand(2/(obs_locy - sind(60)));
I = abs(x(:,1))>=phi;
x(:,1) = x(:,1)-phi;
I = ones(size(I));
y = 0.5*I.*( (mean(x(:,1))*mean(x(:,5)) + (b(1)).*(mean(x(:,1))*mean(x(:,5))-x(:,1).*x(:,3))) .* ...
    (mean(x(:,4)) + (b(2)).*(mean(x(:,4))-x(:,2))) ) + b(3);
end
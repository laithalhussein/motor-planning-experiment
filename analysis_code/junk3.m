close all;

x = normrnd(5,2,1000,1);
y = normrnd(5,2,1000,1);
z = x+y;
mu = 10;
sigma = sqrt(8);
zmin = min(z); zmax = max(z); zrange = range(z);
binw = zrange / 30;
edges = zmin + binw*(0:30);
n = histc(z, edges);
bar(edges, n, 'histc');
hold on
xx = zmin:(zrange/1000):zmax;
plot(xx, binw*length(z)*normpdf(xx,mu,sigma), 'r', 'linewidth', 3);
title('Normal Fit');
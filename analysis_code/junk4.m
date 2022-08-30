
close all;

d1 = normrnd(4,2,[1e3,1]);
d2 = normrnd(0, 3,[1e3,1]);

bin_size = 50;

%get the correct range over which the convolution will rest
z = [d1;d2];
zmin = min(z); 
zmax = max(z); 

bins = linspace(zmin, zmax, (floor(range(z))+1));

[yd1, xd1] = hist(d1, bins);
[yd2, xd2] = hist(d2, bins);

figure; hold on;
histogram(d1,bins);
histogram(d2,bins);

pdf_d1 = yd1/length(d1);
pdf_d2 = yd2/length(d2);
figure; hold on; title('PDFs');
plot(xd2,pdf_d2, 'go');
plot(xd1,pdf_d1, 'ro');

cdf_d1 = cumsum(pdf_d1);
cdf_d2 = cumsum(pdf_d2);
figure; hold on; title('CDFs');
plot(xd1,cdf_d1, 'ro');
plot(xd2,cdf_d2, 'go');

%compute the sum distribution 
pdf_conv = conv(pdf_d1, pdf_d2, 'full');
figure; hold on;
plot(xd1,pdf_d1, 'ro');
plot(xd2,pdf_d2, 'go');

plot(cumsum(pdf_conv));


new_bins = linspace(2*zmin, 2*zmax, length(pdf_conv));
plot(new_bins,pdf_conv, 'k', 'linewidth', 3);

%plot(xd1, yd1, 'o');

function p = boot_test_mean(x1,x2,nboot)
%assumes that samples have the same size

%first compute the mean of the concatinated data
x = [x1(:);x2(:)];
xm = mean(x);

x1_shifted = x1-mean(x1)+xm;
x2_shifted = x2-mean(x2)+xm;

r = length(x1); %number of participants we will resample

x1_boot = nan(nboot,1);
x2_boot = nan(nboot,1);
for k=1:nboot
    observation_boot = randsample(r,r,true);
    x1_boot(k) = mean(x1_shifted(observation_boot));
    x2_boot(k) = mean(x2_shifted(observation_boot));
end

%p value will be the proportion of samples in the null that are greater than the
%observed test statistic
null_dist=x2_boot-x1_boot;
% figure; hold on;
% hist(null_dist,50);
p = sum(null_dist>=mean(x2)-mean(x1))/nboot;
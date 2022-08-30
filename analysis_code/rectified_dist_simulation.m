%In this function we simulate the effects that truncation/rectification have on mean
%and variance of uniform and gaussian distribution

%In both cases, we expect mean to increase, and variance to increase

num_sim = 1000;
u_mean = nan(num_sim, 2);
u_sigma = nan(num_sim, 2);
g_mean = nan(num_sim, 2);
g_sigma = nan(num_sim, 2);

for k=1:num_sim
    %uniform distribution
    u_unrect = unifrnd(-1,1,[1000,1]);
    u_tmp = u_unrect;
    %u_tmp(u_tmp<0)=0;
    u_tmp = abs(u_tmp);
    u_rect = u_tmp;
    
    u_mean(k,1) = mean(u_unrect);
    u_mean(k,2) = mean(u_rect);
    u_sigma(k,1) = std(u_unrect);
    u_sigma(k,2) = std(u_rect);
    
    %gaussian distribution
    g_unrect = normrnd(0,1,[1000,1]);
    g_tmp = g_unrect;
    %g_tmp(g_tmp<0)=0;
    g_tmp = abs(g_tmp);
    g_rect = g_tmp;
    
    g_mean(k,1) = mean(g_unrect);
    g_mean(k,2) = mean(g_rect);
    g_sigma(k,1) = std(g_unrect);
    g_sigma(k,2) = std(g_rect);
    
end

u_ratio = mean(u_sigma(:,2))/mean(u_sigma(:,1))
g_ratio = mean(g_sigma(:,2))/mean(g_sigma(:,1))


function [b_ma, b_po] = determine_model_weights(x,y_ma,y_po,vx,vy,rsub)

%for MA, the naive prediction is to assign equal weights for movement directions asociated with both target (0.5)
%for PO, the naive prediction is to assign equal weight for the movement
%direction that prioritizes obstacle avoidance (the slope of which is the ratio of 2tt:1tt variance) and the movement direction
%that prioritizes movement timing (0 degrees)

%But we should also consider if these weights really should be equal. The data in exp2a indicates that subjects
%actually have a preference in assigning more weight to the obstacle-obstructed target, because the data lie BELOW naive MA prediction
%we can estimate what this weight actually is from the data by calculating the mean of the regression lines that pass through each subject's data point, or a
%regression through all subject data. Then we can adjust the MA prediction
%in exp2a, and use this same regression line in the exp 2b predictions, and since there will be higher weight for the obstacle-obstructed target, then
%the MA prediction in this case will be even further away from the data compared to the 50/50 case

%A minor note is that becuase we have error in both x & y, we could consider a
%deming/PC regression

%the same logic should be applied to the PO predictions. Which is to say,
%regress the data onto [0,...phi,x]

%should we draw CI's on the plots as well? Consider what they would look like...

%To force a regression line to go through a single subject's point: 
%Let's say the subject point is (xi,yi). I will re-center the data with that point as the origin. 
%That is, I subtract xi from every x-value, and yi from every y-value. Now the point is at the origin of the coordinate plane. 
%Then I simply fit a regression line while suppressing the intercept
%Because this is a linear transformation, I can easily back-transform everything afterwards if necessary

[num_subjects, ~] = size(x);

b_ma.sub.dem = nan(1, num_subjects); %holds coefficients for each subject

b_ma.sub.ols = nan(1,num_subjects);
b_ma.sub.ols_inv = nan(1,num_subjects);

% figure; hold on;
% plot(x,y,'.');
for k=1:num_subjects
    
    %first center the data around the response and predictor for that subject
    xtmp = x - x(k);
    ytmp = y_ma - y_ma(k);
    
    %pass them into the regression function
%     [btmp] = eiv_regression(ytmp+mean(ytmp), xtmp+mean(xtmp), vy, vx);
%     b_ma.sub.dem(k) = btmp;
    
    %try regular regression
    b_ma.sub.ols(k) = nlinfit(xtmp,ytmp,@single_dof_ma_model,0.5);
    b_ma.sub.ols_inv(k) = nlinfit(ytmp,xtmp,@single_dof_ma_model,0.5);
    
%     plot(x,x*b_sub(2,k)  + b_sub(1,k),'k'); 
end

% b_ma.avg.dem = eiv_regression(y_ma,x,vy,vx,1);

%try several different intial guesses
ig = [0.5:0.05:0.9];
b_tmp = ig*nan;
for k=1:length(ig)
    b_tmp(k) = nlinfit(x,y_ma,@single_dof_ma_model,ig(k));
end
%if one is different, flag it
if ~all(round(b_tmp,3,'significant')==round(b_tmp(1),3,'significant')), keyboard; end

b_ma.avg.ols = nlinfit(x,y_ma,@single_dof_ma_model,0.5);
b_ma.avg.ols_inv = nlinfit(y_ma,x,@single_dof_ma_model,0.5);

b_ma2 = lsqcurvefit(@single_dof_ma_model,0.5,x,y_ma,0,1);

%% repeat for PO

x_po = [x,rsub];
b_po.avg.ols = nlinfit(x_po,y_po,@single_dof_po_model,[0.3]);

%% plot result
figure; hold on;
plot(x,y_ma,'.');
fit = single_dof_ma_model(b_ma.avg.ols, x);
plot([x],fit,'k');
% plot(x,x*mean(b.sub.ols) + mean(b.sub.ols),'k--');
% plot(x,x*0.5,'color', [1,0,1]/2);

%is the total error negative? Check
err_fit = sum(y_ma - fit);

mse_fit = mean((y_ma-fit).^2);

%make sure the mse of the main fit is lower than adjacent points that could be used for fitting (this will tell us its at a minimum)
fit_ub = single_dof_ma_model(b_ma.avg.ols+0.01, x); %upper bound
fit_lb = single_dof_ma_model(b_ma.avg.ols-0.01, x); %lower bound

mse_ub = mean((y_ma-fit_ub).^2);
mse_lb = mean((y_ma-fit_lb).^2);

%disp(['MSE is ', num2str(mse_fit), ' upper bound is ', num2str(mse_ub), ' and lower bound is ', num2str(mse_lb)]);

%show the mean
%plot(mean(x),mean(y),'p');

%keyboard;

end
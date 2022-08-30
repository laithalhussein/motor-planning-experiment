function [z_nonlin, z_lin] = fit_summary_models_V2(m1, m2, s1, s2, dt_std_gm, st_std_gm, data_title)

%NOTE: m1 should be passed in already flipped!

%in this version we fit the nonlinear and linearized versions of the model:
% (dropping the 0.5), and include offset?
% mu2i = ( K1*(s2i) + (1-K1)*mean(s2) ) * ( K2*mu1i + (1-K2)*mean(mu1) )*( K3*(1/s1i) + (1-K3)*mean(1/s1) )
% =  ( K1*(s2i-mean(s2)) + 1*mean(s2)) * ( K2*(mu1i-mean(mu1)) + 1*mean(mu1) )*( K3*(1/s1i-mean(1/s1) ) +1*mean(1/s1))
% =  [mean(s2)*mean(mu1)*mean(1/s1)]*( K1*(s2i-mean(s2))/mean(s2) + 1)*( K2*(mu1i-mean(mu1))/mean(mu1) + 1)*( K3*(1/s1i-mean(1/s1) )/mean(1/s1) +1)
% = G * (1+A)(1+B)(1+C)  = G* (1+A+B+C+AB+BC+AC+ABC)
% now if A,B,C are all pretty small compared to 1 then:
% ~= G  + (A+B+C)*G

%m1 is vector that hold's each subject's 1-target trial safety margin
%m2 is vector that hold's each subject's 2-target trial IMD
%s is each subject's ratio of SDs

%we analyze the mean-based model, the model individualized for each
%parameter, and then the fully individualized model (5 cases)

%this function should return:
%1) bar graph of SD of prediction error for each case
%2) partial R^2 for each of the three parameters
%3) bar graph of coefficients for K1, K2, and K3

%X = [m1,s];
%X = [m1,s,s*0+sgm];
X = [m1,s1,s2, s1*0+st_std_gm, s2*0+dt_std_gm];
num_subjects = length(m2);

%remove samples that dont satisfy safety margin criteria
% obs_locy = 0.5 * 20 * sind(60);
% phi = atand(2/(obs_locy - sind(60)));
% I = abs(m1)>=phi;
% X(I==0,:) = [];
% m2(I==0,:) = [];

%perhaps bootstrap for K1, K2, and K3?

%make sure to constrain the parameters when fitting
LB = [0,0,0,0,-Inf];
UB = [1,1,1,1,Inf];

%below is the fully individualized model
full_initial_guess = [0.9,0.9,0.9,0.5,0];
z_nonlin.full.params = lsqcurvefit(@full_nonlin_PO_model,full_initial_guess, X, m2, LB, UB);
full_nonlin_model_pred = full_nonlin_PO_model(z_nonlin.full.params, X);
[z_nonlin.full.R2, z_nonlin.full.R2_pval] = calculate_R2(m2,full_nonlin_model_pred,5);
z_nonlin.full.res = m2-full_nonlin_model_pred;
%z_nlin.full.rmse = sqrt(mean((m2-full_nonlin_model_pred).^2));
z_nonlin.full.res_sd = std(z_nonlin.full.res);

%Next we look at the null model (no individualization - everything is based only on means)
null_initial_guess = [0.5,0];
z_nonlin.null.params = lsqcurvefit(@null_nonlin_PO_model,null_initial_guess, X, m2, LB, UB);
null_nonlin_model_pred = null_nonlin_PO_model(z_nonlin.null.params, X);
[z_nonlin.null.R2, z_nonlin.null.R2_pval] = calculate_R2(m2,null_nonlin_model_pred,2);
z_nonlin.null.res = m2-null_nonlin_model_pred;
z_nonlin.null.rmse = sqrt(mean((m2-null_nonlin_model_pred).^2));
%z_nonlin.null.err_ci = std(sqrt((m2-null_nonlin_model_pred).^2))/num_subjects*1.96;
z_nonlin.null.res_sd = std(z_nonlin.null.res);

%Next we look at the s2 individualized model
s2_initial_guess = [0.9,0.5,0];
z_nonlin.s2.params = lsqcurvefit(@s2_nonlin_PO_model,s2_initial_guess, X, m2, LB, UB);
s2_model_pred = s2_nonlin_PO_model(z_nonlin.s2.params, X);
[z_nonlin.s2.R2, z_nonlin.s2.R2_pval] = calculate_R2(m2,s2_model_pred,2);
z_nonlin.s2.res = m2-s2_model_pred;
z_nonlin.s2.rmse = sqrt(mean((m2-s2_model_pred).^2));
z_nonlin.s2.res_sd = std(z_nonlin.s2.res);

%Look at the mu1 individualized model
mu1_initial_guess = [0.9,0.5,0];
z_nonlin.mu1.params = lsqcurvefit(@mu1_nonlin_PO_model,mu1_initial_guess, X, m2, LB, UB);
mu1_model_pred = mu1_nonlin_PO_model(z_nonlin.mu1.params, X);
[z_nonlin.mu1.R2, z_nonlin.mu1.R2_pval] = calculate_R2(m2,mu1_model_pred,2);
z_nonlin.mu1.res = m2-mu1_model_pred;
z_nonlin.mu1.rmse = sqrt(mean((m2-mu1_model_pred).^2));
z_nonlin.mu1.res_sd = std(z_nonlin.mu1.res);

%Look at the s1 individualized model
s1_initial_guess = [0.9,0.5,0];
z_nonlin.s1.params = lsqcurvefit(@s1_nonlin_PO_model,s1_initial_guess, X, m2, LB, UB);
s1_model_pred = s1_nonlin_PO_model(z_nonlin.s1.params, X);
[z_nonlin.s1.R2, z_nonlin.s1.R2_pval] = calculate_R2(m2,s1_model_pred,2);
z_nonlin.s1.res = m2-s1_model_pred;
z_nonlin.s1.rmse = sqrt(mean((m2-s1_model_pred).^2));
z_nonlin.s1.res_sd = std(z_nonlin.s1.res);

%calculate partial R2 of full compared to s2-dropped model
s2_dropped_initial_guess = [0.9,0.9,0.5,-10];
s2_dropped_params = lsqcurvefit(@s2_dropped_nonlin_PO_model,s2_dropped_initial_guess, X, m2, LB, UB);
s2_dropped_model_pred = s2_dropped_nonlin_PO_model(s2_dropped_params, X);
[z_nonlin.s2.pR2, z_nonlin.s2.pR2_pval, z_nonlin.s2.Fstat] = calculate_pR2(m2,s2_dropped_model_pred,full_nonlin_model_pred,...
    length(s2_dropped_initial_guess),length(full_initial_guess));

%calculate partial R2 of full compared to mu1-dropped model
mu1_dropped_initial_guess = [0.9,0.9,0.5,-10];
mu1_dropped_params = lsqcurvefit(@mu1_dropped_nonlin_PO_model,mu1_dropped_initial_guess, X, m2, LB, UB);
mu1_dropped_model_pred = mu1_dropped_nonlin_PO_model(mu1_dropped_params, X);
[z_nonlin.mu1.pR2, z_nonlin.mu1.pR2_pval, z_nonlin.mu1.Fstat] = calculate_pR2(m2,mu1_dropped_model_pred,full_nonlin_model_pred,...
    length(mu1_dropped_initial_guess),length(full_initial_guess));

%calculate partial R2 of full compared to s1-dropped model
s1_dropped_initial_guess = [0.9,0.9,0.5,-10];
s1_dropped_params = lsqcurvefit(@s1_dropped_nonlin_PO_model,s1_dropped_initial_guess, X, m2, LB, UB);
s1_dropped_model_pred = s1_dropped_nonlin_PO_model(s1_dropped_params, X);
[z_nonlin.s1.pR2, z_nonlin.s1.pR2_pval, z_nonlin.s1.Fstat] = calculate_pR2(m2,s1_dropped_model_pred,full_nonlin_model_pred,...
    length(s1_dropped_initial_guess),length(full_initial_guess));

%% plot SD of the residuals
figure; hold on;
ylabel('SD of individual participant error (deg)');

xbar_offset = 5;
green = [0 0.5 0];

%if the S1 SD and null SD are about the same, decrease S1 for visual clarity

sd_res_all.nonlin = [z_nonlin.null.res_sd; z_nonlin.mu1.res_sd; z_nonlin.s1.res_sd; z_nonlin.s2.res_sd; z_nonlin.full.res_sd];

bar([1,2,3,4,5],sd_res_all.nonlin,0.25,'FaceColor', green, 'EdgeColor', 'None'); hold on;

xticks([1,2,3,4,5,xbar_offset+[4:5]]);
xticklabels({'Group mean PO model','Individuation of M1','Individuation of S1','Individuation of S2', 'Fully Individuated PO model'});
xtickangle(45)

ylim([0,4.5]);
set(gca, 'YTick', [0:0.5:4],'YTickLabel', {'0', '', '1', '', '2', '', '3', '', '4', ''});
htmp=gca; 
htmp.XAxis.TickLength = [0, 0];

title(['residual SD in ',data_title]);

%% plot partial R^2

figure; hold on;
ylabel('partial R^2 for PO model parameters');

pR2_all.nonlin = [z_nonlin.mu1.pR2; z_nonlin.s1.pR2; z_nonlin.s2.pR2];

bar([1,2,3],pR2_all.nonlin,0.5,'FaceColor', 'k', 'EdgeColor', 'None');

xticks([1,2,3]);
xticklabels({'\mu_1','\sigma_1','\sigma_2'});

ylim([0,1]);
set(gca, 'YTick', [0:0.5:1],'YTickLabel', {'0', '0.5', '1'});
htmp=gca; 
htmp.XAxis.TickLength = [0, 0];

title(['partial R^2 for ',data_title]);


%% repeat for linear model

z_lin=NaN;


%%
% figure; hold on;
% title('Exp 2b model errors');
% 
% hb1 = barwitherr( [z_nonlin.null.err_ci, z_nonlin.alpha.err_ci, z_nonlin.full.err_ci, z_nonlin.beta.err_ci],...
%     [z_nonlin.null.rmse, z_nonlin.alpha.rmse, z_nonlin.full.rmse, z_nonlin.beta.rmse] );
% 
% set(gca, 'XTick', [1,2,3,4],'XTickLabel', {'Null model', 'Mean model', 'Full model', 'Variability model'});
% ylabel('RMSE (deg)');




% (1) no individuation (as a baseline)
% (2) all-measurement individuation (using personalized values of all three variables: u1, s1, & s2
% (3) single-measurement individuation (using personalized values one variable at a time)
%            (3a) single-measurement individuation for u1 (use group-mean values for s1 & s2)
%            (3b) single-measurement individuation for s1 (use group-mean values for u1 & s2)
%            (3c) single-measurement individuation for s2 (use group-mean values for u1 & s1)

%the raw null model is the no individuation case (means only)
%the raw full model is using full individuation
%the raw alpha model allows individuation for u1 only
%the tol beta model does individuation for s2


end
function [z_raw, z_ref, z_tol, z_indiv] = fit_summary_models(m1, m2, s1, s2, dt_std_gm, st_std_gm)

%m1 is vector that hold's each subject's 1-target trial safety margin
%m2 is vector that hold's each subject's 2-target trial IMD
%s is each subject's ratio of SDs

%fit the following models for a summary analysis:
%null model: m2_i = mean(m1)*mean(s) + offset
%full (individualized) model: m2_i = ( mean(m1) +alpha*(m1_i-mean(m1) ) * (mean(s) + beta*(s_i-mean(s) ) + offset

s = s2./s1;
%X = [m1,s];
sgm = dt_std_gm/st_std_gm;
X = [m1,s,s*0+sgm];
num_subjects = length(m2);

%remove samples that dont satisfy safety margin criteria
% obs_locy = 0.5 * 20 * sind(60);
% phi = atand(2/(obs_locy - sind(60)));
% I = abs(m1)>=phi;
% X(I==0,:) = [];
% m2(I==0,:) = [];

z_raw.full.params = nlinfit(X,m2,@full_PO_model,[0.9,0.9,30]);
full_model_pred = full_PO_model(z_raw.full.params, X);
[z_raw.full.R2, z_raw.full.R2_pval] = calculate_R2(m2,full_model_pred,3);
z_raw.full.res = m2-full_model_pred;
z_raw.full.rmse = sqrt(mean((m2-full_model_pred).^2));
z_raw.full.err_ci = std(sqrt((m2-full_model_pred).^2))/num_subjects*1.96;

z_indiv.full.params = nlinfit(X,m2,@indiv_PO_model,[10]);
indiv_model_pred = indiv_PO_model(z_indiv.full.params, X);
[z_indiv.full.R2, z_indiv.full.R2_pval] = calculate_R2(m2,indiv_model_pred,3);
z_indiv.full.res = m2-indiv_model_pred;
z_indiv.full.rmse = sqrt(mean((m2-indiv_model_pred).^2));
z_indiv.full.err_ci = std(sqrt((m2-indiv_model_pred).^2))/num_subjects*1.96;

z_raw.null.params = nlinfit(X,m2,@null_PO_model,30);
null_model_pred = null_PO_model(z_raw.null.params, X);
[z_raw.null.R2, z_raw.null.R2_pval] = calculate_R2(m2,null_model_pred,1);
z_raw.null.res = m2-null_model_pred;
z_raw.null.rmse = sqrt(mean((m2-null_model_pred).^2));
z_raw.null.err_ci = std(sqrt((m2-null_model_pred).^2))/num_subjects*1.96;

z_raw.alpha.params = nlinfit(X,m2,@alpha_PO_model,[0.5,30]);
alpha_model_pred = alpha_PO_model(z_raw.alpha.params, X);
[z_raw.alpha.R2, z_raw.alpha.R2_pval] = calculate_R2(m2,alpha_model_pred,2);
z_raw.alpha.res = m2-alpha_model_pred;
z_raw.alpha.rmse = sqrt(mean((m2-alpha_model_pred).^2));
z_raw.alpha.err_ci = std(sqrt((m2-alpha_model_pred).^2))/num_subjects*1.96;

z_raw.beta.params = nlinfit(X,m2,@beta_PO_model,[0.5,30]);
beta_model_pred = beta_PO_model(z_raw.beta.params, X);
[z_raw.beta.R2, z_raw.beta.R2_pval] = calculate_R2(m2,beta_model_pred,2);
z_raw.beta.res = m2-beta_model_pred;
z_raw.beta.rmse = sqrt(mean((m2-beta_model_pred).^2));
z_raw.beta.err_ci = std(sqrt((m2-beta_model_pred).^2))/num_subjects*1.96;

z_raw.constrained.params = nlinfit(X,m2,@constrained_PO_model,[30]);
constrained_model_pred = constrained_PO_model(z_raw.constrained.params, X);
[z_raw.constrained.R2, z_raw.constrained.R2_pval] = calculate_R2(m2,constrained_model_pred,1);
z_raw.constrained.res = m2-full_model_pred;
z_raw.constrained.rmse = sqrt(mean((m2-constrained_model_pred).^2));
z_raw.constrained.err_ci = std(sqrt((m2-constrained_model_pred).^2))/num_subjects*1.96;

%calculate partial R2 of full compared to alpha only model
[z_raw.alpha.pR2, z_raw.alpha.pR2_pval, z_raw.alpha.Fstat] = calculate_pR2(m2,alpha_model_pred,full_model_pred,2,3);

%calculate partial R2 of full compared to beta only model
[z_raw.beta.pR2, z_raw.beta.pR2_pval, z_raw.beta.Fstat] = calculate_pR2(m2,beta_model_pred,full_model_pred,2,3);

%plot RMSE
figure; hold on;
title('Exp 2b model errors');

hb1 = barwitherr( [z_raw.null.err_ci, z_raw.alpha.err_ci, z_raw.full.err_ci, z_raw.beta.err_ci],...
    [z_raw.null.rmse, z_raw.alpha.rmse, z_raw.full.rmse, z_raw.beta.rmse] );

set(gca, 'XTick', [1,2,3,4],'XTickLabel', {'Null model', 'Mean model', 'Full model', 'Variability model'});
ylabel('RMSE (deg)');

%% repeat the analysis for the refined PO model

%full refined model
z_ref.full.params = nlinfit(X,m2,@full_ref_PO_model,[1,0.9,0.9,30]);
full_ref_model_pred = full_ref_PO_model(z_ref.full.params, X);
[z_ref.full.R2, z_ref.full.R2_pval] = calculate_R2(m2,full_ref_model_pred,4);
z_ref.full.res = m2-full_ref_model_pred;
z_ref.full.rmse = sqrt(mean((m2-full_ref_model_pred).^2));
z_ref.full.err_ci = std(sqrt((m2-full_ref_model_pred).^2))/num_subjects*1.96;

%null refined model
z_ref.null.params = nlinfit(X,m2,@null_ref_PO_model,[1,30]);
null_ref_model_pred = null_ref_PO_model(z_ref.null.params, X);
[z_ref.null.R2, z_ref.null.R2_pval] = calculate_R2(m2,null_ref_model_pred,2);
z_ref.null.res = m2-null_ref_model_pred;
z_ref.null.rmse = sqrt(mean((m2-null_ref_model_pred).^2));
z_ref.null.err_ci = std(sqrt((m2-null_ref_model_pred).^2))/num_subjects*1.96;

%alpha-only refined model
z_ref.alpha.params = nlinfit(X,m2,@alpha_ref_PO_model,[1,0.5,30]);
alpha_ref_model_pred = alpha_ref_PO_model(z_ref.alpha.params, X);
[z_ref.alpha.R2, z_ref.alpha.R2_pval] = calculate_R2(m2,alpha_ref_model_pred,3);
z_ref.alpha.res = m2-alpha_ref_model_pred;
z_ref.alpha.rmse = sqrt(mean((m2-alpha_ref_model_pred).^2));
z_ref.alpha.err_ci = std(sqrt((m2-alpha_ref_model_pred).^2))/num_subjects*1.96;

%beta-only refined model
z_ref.beta.params = nlinfit(X,m2,@beta_ref_PO_model,[1,0.5,30]);
beta_ref_model_pred = beta_ref_PO_model(z_ref.beta.params, X);
[z_ref.beta.R2, z_ref.beta.R2_pval] = calculate_R2(m2,beta_ref_model_pred,3);
z_ref.beta.res = m2-beta_ref_model_pred;
z_ref.beta.rmse = sqrt(mean((m2-beta_ref_model_pred).^2));
z_ref.beta.err_ci = std(sqrt((m2-beta_ref_model_pred).^2))/num_subjects*1.96;

%calculate partial R2 of full refined model compared to alpha-only refined model
[z_ref.alpha.pR2, z_ref.alpha.pR2_pval, z_ref.alpha.Fstat] = calculate_pR2(m2,alpha_ref_model_pred,full_ref_model_pred,3,4);

%calculate partial R2 of full refined model compared to beta-only refined model
[z_ref.beta.pR2, z_ref.beta.pR2_pval, z_ref.beta.Fstat] = calculate_pR2(m2,beta_ref_model_pred,full_ref_model_pred,3,4);


%% fit data onto risk tolerance model:
%m2i = 0.5*( m1/s1 + alpha*( m1/s1 - m1i/s1i ) * ( s2 + beta*(s2 - s2i) ) + offset

X2 = [m1,s1,s2, s1*0+st_std_gm, s2*0+dt_std_gm];
% X2(I==0,:) = [];

%fit the full variant
z_tol.full.params = nlinfit(X2,m2,@full_tol_PO_model,[0.5,0.5,15]);
full_tol_model_pred = full_tol_PO_model(z_tol.full.params, X2);
[z_tol.full.R2, z_tol.full.R2_pval] = calculate_R2(m2,full_tol_model_pred,3);
z_tol.full.res = m2-full_tol_model_pred;
z_tol.full.rmse = sqrt(mean((m2-full_tol_model_pred).^2));
z_tol.full.err_ci = std(sqrt((m2-full_tol_model_pred).^2))/num_subjects*1.96;

%fit to null variant
z_tol.null.params = nlinfit(X2,m2,@null_tol_PO_model,[15]);
null_tol_model_pred = null_tol_PO_model(z_tol.null.params, X2);
[z_tol.null.R2, z_tol.null.R2_pval] = calculate_R2(m2,null_tol_model_pred,1);
z_tol.null.res = m2-null_tol_model_pred;
z_tol.null.rmse = sqrt(mean((m2-null_tol_model_pred).^2));
z_tol.null.err_ci = std(sqrt((m2-null_tol_model_pred).^2))/num_subjects*1.96;

%fit to alpha-only variant
z_tol.alpha.params = nlinfit(X2,m2,@alpha_tol_PO_model,[0.5,15]);
alpha_tol_model_pred = alpha_tol_PO_model(z_tol.alpha.params, X2);
[z_tol.alpha.R2, z_tol.alpha.R2_pval] = calculate_R2(m2,alpha_tol_model_pred,2);
z_tol.alpha.res = m2-alpha_tol_model_pred;
z_tol.alpha.rmse = sqrt(mean((m2-alpha_tol_model_pred).^2));
z_tol.alpha.err_ci = std(sqrt((m2-alpha_tol_model_pred).^2))/num_subjects*1.96;

%fit to beta-only variant
z_tol.beta.params = nlinfit(X2,m2,@beta_tol_PO_model,[0.5,15]);
beta_tol_model_pred = beta_tol_PO_model(z_tol.beta.params, X2);
[z_tol.beta.R2, z_tol.beta.R2_pval] = calculate_R2(m2,beta_tol_model_pred,2);
z_tol.beta.res = m2-beta_tol_model_pred;
z_tol.beta.rmse = sqrt(mean((m2-beta_tol_model_pred).^2));
z_tol.beta.err_ci = std(sqrt((m2-beta_tol_model_pred).^2))/num_subjects*1.96;

%calculate partial R2 of full tolerance model compared to alpha-only refined model
[z_tol.alpha.pR2, z_tol.alpha.pR2_pval, z_tol.alpha.Fstat] = calculate_pR2(m2,alpha_tol_model_pred,full_tol_model_pred,2,3);

%calculate partial R2 of full tolerance model compared to beta-only refined model
[z_tol.beta.pR2, z_tol.beta.pR2_pval, z_tol.beta.Fstat] = calculate_pR2(m2,beta_tol_model_pred,full_tol_model_pred,2,3);

% figure; hold on;
% title('PO model variant errors');
% 
% hb1 = barwitherr( [z_raw.full.err_ci, z_tol.full.err_ci],...
%     [z_raw.full.rmse, z_tol.full.rmse] );
% 
% set(gca, 'XTick', [1,2],'XTickLabel', {'Full model', 'Risk tolerance model'});
% ylabel('RMSE (deg)');


% (1) no individuation (as a baseline)
% (2) all-measurement individuation (using personalized values of all three variables: u1, s1, & s2
% (3) single-measurement individuation (using personalized values one variable at a time)
%            (3a) single-measurement individuation for u1 (use group-mean values for s1 & s2)
%            (3b) single-measurement individuation for s1 (use group-mean values for u1 & s2)
%            (3c) single-measurement individuation for s2 (use group-mean values for u1 & s1)

figure; hold on;
title('Exp 2b model errors');

%the raw null model is the no individuation case (means only)
%the raw full model is using full individuation
%the raw alpha model allows individuation for u1 only
%the tol beta model does individuation for s2


hb1 = barwitherr( [z_raw.null.err_ci, z_raw.alpha.err_ci, z_raw.full.err_ci, z_raw.beta.err_ci],...
    [z_raw.null.rmse, z_raw.alpha.rmse, z_raw.full.rmse, z_raw.beta.rmse] );

set(gca, 'XTick', [1,2,3,4],'XTickLabel', {'No individuation', 'Full individuation', 'Indiviudation for u1',...
    'Indiviudation for s1', 'Indiviudation for s2'});
ylabel('RMSE (deg)');


end
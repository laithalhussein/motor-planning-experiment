close all;

%make figures in response to R2

%% first, extract the MA and PO predictions for when the safety margin size is doubled and halved

sm_changes = [0.5,1,2]; %look at it when its the same, or doubled, or halved
sm_fields = {'sm1','sm2','sm3'};
close all;
clear ma po po_scaled;
%sm, phi, obs_ratio_all, obs_ratio_gm, slope
for k=1:length(sm_changes)
    cfld = sm_fields{k};
    po.(cfld) = determine_po_pred(sm_diff*sm_changes(k), phiA, obs_ratio_sub, obs_ratio_gm, 0.5);
    %po_scaled.(cfld) = determine_po_pred_scaled(sm_diff*sm_changes(k), phiA, obs_ratio_sub, obs_ratio_gm, 0.5);
    ma.(cfld) = determine_ma_pred(sm_diff*sm_changes(k), obs_nb_bias.comb, 0.5);
end
close all;
figure; hold on;

plot([1], dtt_obs_avg.expt2b, 's', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'markersize', 8);
display_xy_error_V2([1],dtt_obs_avg.expt2b, [], dtt_obs_err.expt2b, 'k');

hline(mean(po.sm1.gm), 'linestyle', '--', 'color', 'g');
hline(mean(ma.sm1.sub), 'linestyle', '--', 'color', [0.5,0,0.5]);

hline(mean(po.sm2.gm), 'linestyle', '-', 'color', 'g');
hline(mean(ma.sm2.sub), 'linestyle', '-', 'color', [0.5,0,0.5]);

hline(mean(po.sm3.gm), 'linestyle', '-.', 'color', 'g');
hline(mean(ma.sm3.sub), 'linestyle', '-.', 'color', [0.5,0,0.5]);

title('Expt2b: Effect of obstacle-present 1-target trial safety margin size');
set(gca, 'XTick', [1],'XTickLabel', { 'obstacle-present 2-target trial data'});
ylabel('IMD (deg)');
ylim([-22,22]);
set(gca, 'YTick', [-20:5:20],'YTickLabel', {'-20', '', '-10', '', '0', '', '10', '', '20'});

%% repeat the above, but this time, on the x-axis, plot the percent change from population-averaged safety margin on x-axis and
%on y-axis, plot the expected MA and PO prediction

sm_changes = [0.5:0.01:2]; %min and max is 58% to ~175%
ma_smc = nan(length(sm_changes),1);
po_smc = nan(length(sm_changes),1);
ma_smc_err = nan(length(sm_changes),1);
po_smc_err = nan(length(sm_changes),1);

close all;
clear ma po;
for k=1:length(sm_changes)
    po_tmp = determine_po_pred(sm_diff*sm_changes(k), phiA, obs_ratio_sub, obs_ratio_gm, 1);
    ma_tmp = determine_ma_pred(sm_diff*sm_changes(k), obs_nb_bias.comb);
    po_smc(k) = mean(po_tmp.sub);
    po_smc_err(k) = std(po_tmp.sub)/sqrt(length(po_tmp.sub))*1.96;
    ma_smc(k) = mean(ma_tmp.sub);
    ma_smc_err(k) = std(ma_tmp.sub)/sqrt(length(ma_tmp.sub))*1.96;
end
close all;
figure; hold on;

xlabel('Percentage of population-averaged 1TT safety margin');
ylabel('Predicted IMD on obstacle-present 2-target trials (deg)');

plot(sm_changes*100, po_smc, 'color', 'g', 'linewidth', 2);
display_xy_error_V2(sm_changes*100, po_smc, [], po_smc_err, 'g');
plot(sm_changes*100, ma_smc, 'color', [0.5,0,0.5], 'linewidth', 2);
display_xy_error_V2(sm_changes*100, ma_smc, [], ma_smc_err, [0.5,0,0.5]);

ylim([-22,22]);

%% Lets go further and produce plots similar to Figs 5b/e and show the inidividual data, but with the lines of the predictions based
%on halved and doubled safety margins

sm_changes = [0.5,1,2]; %look at it when its the same, or doubled, or halved
close all;

figure; hold on;
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    display_PO_diff_flipped_R2_sm(sm_diff,dt_diff, 1, phiA, [], 0.66, 1-0.722, obs_nb_bias,...
    obs_ratio_sub, obs_ratio_gm, 0, sm_changes);

%% repeat all the above analyses, but this time, manipulate the ratios instead
ratio_changes = [0.5,1,2];

sm_fields = {'rat1','rat2','rat3'};
close all;
clear ma po;
for k=1:length(sm_changes)
    cfld = sm_fields{k};
    po.(cfld) = determine_po_pred(sm_diff, phiA, obs_ratio_sub*ratio_changes(k), obs_ratio_gm*ratio_changes(k), 1);
    ma.(cfld) = determine_ma_pred(sm_diff, obs_nb_bias.comb);
end
close all;
figure; hold on;

plot([1], dtt_obs_avg.expt2b, 's', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'markersize', 8);
display_xy_error_V2([1],dtt_obs_avg.expt2b, [], dtt_obs_err.expt2b, 'k');

hline(mean(po.rat1.sub), 'linestyle', '--', 'color', 'g');
hline(mean(ma.rat1.sub), 'linestyle', '--', 'color', [0.5,0,0.5]);

hline(mean(po.rat2.sub), 'linestyle', '-', 'color', 'g');
hline(mean(ma.rat2.sub), 'linestyle', '-', 'color', [0.5,0,0.5]);

hline(mean(po.rat3.sub), 'linestyle', '-.', 'color', 'g');
hline(mean(ma.rat3.sub), 'linestyle', '-.', 'color', [0.5,0,0.5]);

title('Expt2b: Effect of Eq. 2 slope');
set(gca, 'XTick', [1],'XTickLabel', { 'obstacle-present 2-target trial data'});
ylabel('IMD (deg)');
ylim([-22,22]);
set(gca, 'YTick', [-20:5:20],'YTickLabel', {'-20', '', '-10', '', '0', '', '10', '', '20'});

%% repeat the above, but this time, on the x-axis, plot the percent change from population-averaged variability ratio on x-axis and
%on y-axis, plot the expected MA and PO prediction

rat_changes = [0.5:0.01:2]; %min and max is 58% to ~175%
ma_rc = nan(length(rat_changes),1);
po_rc = nan(length(rat_changes),1);
ma_rc_err = nan(length(rat_changes),1);
po_rc_err = nan(length(rat_changes),1);

close all;
clear ma po;
for k=1:length(rat_changes)
    po_tmp = determine_po_pred(sm_diff, phiA, obs_ratio_sub*rat_changes(k), obs_ratio_gm*rat_changes(k), 1);
    ma_tmp = determine_ma_pred(sm_diff, obs_nb_bias.comb);
    po_rc(k) = mean(po_tmp.sub);
    po_rc_err(k) = std(po_tmp.sub)/sqrt(length(po_tmp.sub))*1.96;
    ma_rc(k) = mean(ma_tmp.sub);
    ma_rc_err(k) = std(ma_tmp.sub)/sqrt(length(ma_tmp.sub))*1.96;
end
close all;
figure; hold on;

xlabel('Percentage of population-averaged variability ratio');
ylabel('Predicted IMD on obstacle-present 2-target trials (deg)');

plot(rat_changes*100, po_rc, 'color', 'g', 'linewidth', 2);
display_xy_error_V2(rat_changes*100, po_rc, [], po_rc_err, 'g');
plot(rat_changes*100, ma_rc, 'color', [0.5,0,0.5], 'linewidth', 2);
display_xy_error_V2(rat_changes*100, ma_rc, [], ma_rc_err, [0.5,0,0.5]);

ylim([-20,20]);

%% Lets go further and produce plots similar to Figs 5b/e and show the inidividual data, but with the lines of the predictions based
%on halved and doubled safety margins

rat_changes = [0.5,1,2]; %look at it when its the same, or doubled, or halved
close all;

%NOTE: Not sure if I should actually show this
figure; hold on;
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    display_PO_diff_flipped_R2_rat(sm_diff,dt_diff, 1, phiA, [], 0.66, 1-0.722, obs_nb_bias,...
    obs_ratio_sub, obs_ratio_gm, 0, rat_changes);

%% what do you get when you fit the refined slopes on Expt 2b first
%in the same figure, plot the 33% and 66% weighting scheme
figure; hold on;
title('exp 2b predictions');
%the PO weight should be 1-0.722
[~, pmodel, smodel, r2_model, corr_pred, p_corr_pred, dt_diff_ma_shifted, dt_diff_po_shifted, po_pred, ma_pred] = ...
    display_PO_diff_flipped_0623_2020a_vref( sm_diff, dt_diff, 1, phiA, [], 0.685, 0.22, obs_nb_bias, obs_ratio_sub, obs_ratio_gm,...
    ['Expt 2b: RMSEs based on fits']); %Note that this returns the RMSEs only for the fits based on refined coefficient

%% use the Expt 2b fits to determine the refined RMSEs for Expt 2a
figure; hold on; title('Expt 2a: RMSEs based on predictions');
%the PO weight should be 1-0.722
[~, pmodel, ~, ~, ~, ~, ~, ~, exp2a.po.(cfld), exp2a.ma.(cfld)] = ...
    display_PO_diff_exp2a_vref( exp2a.sm_diff, exp2a.dt_diff, -1, phiA, [], 0.685, 0.22, exp2a.obs_nb_bias, exp2a.std_obs_ratio,...
    exp2a.obs_ratio_gm, {'Expt 2a: RMSEs based on predictions'});
hold off;

%% check 33% vs 66% weighting and then 66% vs 33% weighting

%the MA weight is for the obstacle obstructed target
%the PO weight is for the obstacle avoidance

%first try 66% for obstacle-obstructed target, and 66% for OA
figure; hold on; title('66% weighting');
[~, pmodel66, ~, ~, ~, ~, ~, ~, po.(cfld), ma.(cfld)] = ...
    display_PO_diff_flipped_0623_2020a_vref( sm_diff, dt_diff, 1, phiA, [], 0.66, 0.66, obs_nb_bias, obs_ratio_sub, obs_ratio_gm, {'Expt 2b: 66% weighting'});
hold off;

%try 33% for obstacle-obstructed target, and 33% for OA
figure; hold on; title('33% weighting');
[~, pmodel33, ~, ~, ~, ~, ~, ~, po.(cfld), ma.(cfld)] = ...
    display_PO_diff_flipped_0623_2020a_vref( sm_diff, dt_diff, 1, phiA, [], 0.33, 0.33, obs_nb_bias, obs_ratio_sub, obs_ratio_gm, {'Expt 2b: 33% weighting'});
hold off;
%% check MA predictions in Expt 2a

%look at original (combined data for each subject)

%safety margins
figure; hold on;
plot([1:exp2a.num_subjects], exp2a.sm_diff, 'ko', 'markersize', 8, 'displayname', 'safety margin');

%recover the original movement direction
plot([1:exp2a.num_subjects], exp2a.sm_diff-phiA, 'rx', 'markersize', 8, 'displayname', 'movement direction');

%get each subject's MA prediction with the safety margin
sm_ma_pred = determine_ma_pred_expt2a_V2(exp2a.sm_diff, exp2a.obs_nb_bias.comb, phiA); %this is incorrect

%get each subject's MA prediction with the original 1TT IMD
md_ma_pred = determine_ma_pred_expt2a_V2(exp2a.sm_diff, exp2a.obs_nb_bias.comb, phiA);
hline(mean(md_ma_pred.sub), 'linestyle', '-', 'color', purple, 'displayname', 'MA prediction with unchanged safety margin');

%get each subject's MA prediction with a halved original 1TT IMD
md2_ma_pred = determine_ma_pred_expt2a_V2((exp2a.sm_diff)*0.5, exp2a.obs_nb_bias.comb, phiA);
hline(mean(md2_ma_pred.sub), 'linestyle', '--', 'color', purple, 'displayname', 'MA prediction with halved safety margin');

legend('show');

title('1-target trial data in Expt 2a');
xlabel('Subject');
ylabel('deg');

%% look at effect of safety margin for Expt 2a
sm_changes = [0.5,1,2]; %look at it when its the same, or doubled, or halved
sm_fields = {'sm1','sm2','sm3'};
close all;
clear ma po po_scaled;
for k=1:length(sm_changes)
    cfld = sm_fields{k};
    po.(cfld) = determine_po_pred_expt2a(exp2a.sm_diff*sm_changes(k), phiA, exp2a.std_obs_ratio, exp2a.obs_ratio_gm, 0.5);
    po_scaled.(cfld) = determine_po_pred_expt2a_scaled(exp2a.sm_diff*sm_changes(k), phiA, exp2a.std_obs_ratio, exp2a.obs_ratio_gm, 0.5);
    ma.(cfld) = determine_ma_pred_expt2a_V2(exp2a.sm_diff*sm_changes(k), exp2a.obs_nb_bias.comb, phiA, 0.5);
end
close all;
figure; hold on;

plot([1], dtt_obs_avg.expt2a, 's', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'markersize', 8);
display_xy_error_V2([1],dtt_obs_avg.expt2a, [], dtt_obs_err.expt2a, 'k');

hline(mean(po.sm1.gm), 'linestyle', '--', 'color', 'g');
hline(mean(ma.sm1.sub), 'linestyle', '--', 'color', [0.5,0,0.5]);

hline(mean(po.sm2.gm), 'linestyle', '-', 'color', 'g');
hline(mean(ma.sm2.sub), 'linestyle', '-', 'color', [0.5,0,0.5]);

hline(mean(po.sm3.gm), 'linestyle', '-.', 'color', 'g');
hline(mean(ma.sm3.sub), 'linestyle', '-.', 'color', [0.5,0,0.5]);

title('Expt2a: Effect of obstacle-present 1-target trial safety margin size');
set(gca, 'XTick', [1],'XTickLabel', { 'obstacle-present 2-target trial data'});
ylabel('IMD (deg)');
ylim([-5,32]);
set(gca, 'YTick', [0:5:30],'YTickLabel', {'0', '', '10', '', '20', '', '30'});

%% perform percent change of safety margin analysis for Expt 2a

sm_changes = [0.5:0.01:2]; %min and max is ?? to ??
ma_smc = nan(length(sm_changes),1);
po_smc = nan(length(sm_changes),1);
ma_smc_err = nan(length(sm_changes),1);
po_smc_err = nan(length(sm_changes),1);

close all;
for k=1:length(sm_changes)
    po_tmp = determine_po_pred_expt2a(exp2a.sm_diff*sm_changes(k), phiA, exp2a.std_obs_ratio, exp2a.obs_ratio_gm, -1);
    ma_tmp = determine_ma_pred_expt2a_V2(exp2a.sm_diff*sm_changes(k), exp2a.obs_nb_bias.comb, phiA);
    po_smc(k) = mean(po_tmp.sub);
    po_smc_err(k) = std(po_tmp.sub)/sqrt(length(po_tmp.sub))*1.96;
    ma_smc(k) = mean(ma_tmp.sub);
    ma_smc_err(k) = std(ma_tmp.sub)/sqrt(length(ma_tmp.sub))*1.96;
end
close all;
figure; hold on;

xlabel('Percentage of population-averaged 1TT safety margin');
ylabel('Predicted IMD on obstacle-present 2-target trials (deg)');

plot(sm_changes*100, po_smc, 'color', 'g', 'linewidth', 2);
h1=display_xy_error_V2(sm_changes*100, po_smc, [], po_smc_err, 'g');
plot(sm_changes*100, ma_smc, 'color', [0.5,0,0.5], 'linewidth', 2);
h2=display_xy_error_V2(sm_changes*100, ma_smc, [], ma_smc_err, [0.5,0,0.5]);

ylim([-5,35]);

%% Lets go further and produce plots similar to Figs 5b/e and show the inidividual data, but with the lines of the predictions based
%on halved and doubled safety margins

sm_changes = [0.5,1,2]; %look at it when its the same, or doubled, or halved
close all;

figure; hold on;
[~, ~, ~, ~, ~, ~, ~, ~, exp2a.po.(cfld), exp2a.ma.(cfld)] = ...
    display_PO_diff_exp2a_R2_sm( exp2a.sm_diff, exp2a.dt_diff, -1, phiA, [], 0.64, 0.22, exp2a.obs_nb_bias,...
    exp2a.std_obs_ratio, exp2a.obs_ratio_gm, sm_changes);

%% manipulate the ratios instead for Expt 2a
ratio_changes = [0.5,1,2];

rat_fields = {'rat1','rat2','rat3'};
close all;
clear ma po;
for k=1:length(sm_changes)
    cfld = rat_fields{k};
    po.(cfld) = determine_po_pred_expt2a(exp2a.sm_diff, phiA, exp2a.std_obs_ratio*ratio_changes(k), exp2a.obs_ratio_gm*ratio_changes(k), -1);
    ma.(cfld) = determine_ma_pred_expt2a_V2(exp2a.sm_diff, exp2a.obs_nb_bias.comb, phiA);
end
close all;
figure; hold on;

plot([1], dtt_obs_avg.expt2a, 's', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'markersize', 8);
display_xy_error_V2([1],dtt_obs_avg.expt2a, [], dtt_obs_err.expt2a, 'k');

hline(mean(po.rat1.sub), 'linestyle', '--', 'color', 'g');
hline(mean(ma.rat1.sub), 'linestyle', '--', 'color', [0.5,0,0.5]);

hline(mean(po.rat2.sub), 'linestyle', '-', 'color', 'g');
hline(mean(ma.rat2.sub), 'linestyle', '-', 'color', [0.5,0,0.5]);

hline(mean(po.rat3.sub), 'linestyle', '-.', 'color', 'g');
hline(mean(ma.rat3.sub), 'linestyle', '-.', 'color', [0.5,0,0.5]);

title('Expt2a: Effect of Eq. 2 slope');
set(gca, 'XTick', [1],'XTickLabel', { 'obstacle-present 2-target trial data'});
ylabel('IMD (deg)');
%ylim([-0,15]);
ylim([-5,30]);
%set(gca, 'YTick', [-5:5:20],'YTickLabel', {'-20', '', '-10', '', '0', '', '10', '', '20'});

%% repeat the above, but this time, on the x-axis, plot the percent change from population-averaged variability ratio on x-axis and
%on y-axis, plot the expected MA and PO prediction

rat_changes = [0.5:0.01:2]; %min and max is 58% to ~175%
ma_rc = nan(length(rat_changes),1);
po_rc = nan(length(rat_changes),1);
ma_rc_err = nan(length(rat_changes),1);
po_rc_err = nan(length(rat_changes),1);

close all;
clear ma po;
for k=1:length(rat_changes)
    po_tmp = determine_po_pred_expt2a(exp2a.sm_diff, phiA, exp2a.std_obs_ratio*rat_changes(k), exp2a.obs_ratio_gm*rat_changes(k), -1);
    ma_tmp = determine_ma_pred_expt2a_V2(exp2a.sm_diff, exp2a.obs_nb_bias.comb, phiA);
    po_rc(k) = mean(po_tmp.sub);
    po_rc_err(k) = std(po_tmp.sub)/sqrt(length(po_tmp.sub))*1.96;
    ma_rc(k) = mean(ma_tmp.sub);
    ma_rc_err(k) = std(ma_tmp.sub)/sqrt(length(ma_tmp.sub))*1.96;
end
close all;
figure; hold on;

xlabel('Percentage of population-averaged variability ratio');
ylabel('Predicted IMD on obstacle-present 2-target trials (deg)');

plot(rat_changes*100, po_rc, 'color', 'g', 'linewidth', 2);
display_xy_error_V2(rat_changes*100, po_rc, [], po_rc_err, 'g');
ylim([-20,25]);

figure; hold on;
plot(rat_changes*100, ma_rc, 'color', [0.5,0,0.5], 'linewidth', 2);
display_xy_error_V2(rat_changes*100, ma_rc, [], ma_rc_err, [0.5,0,0.5]);
ylim([-20,25]);

%% Show the individual data plot for changes in var ratio in Expt 2a

rat_changes = [0.5,1,2]; %look at it when its the same, or doubled, or halved
close all;

%NOTE: Not sure if I should actually show this
figure; hold on;
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    display_PO_diff_exp2a_R2_rat( exp2a.sm_diff, exp2a.dt_diff, -1, phiA, [], 0.64, 0.22, exp2a.obs_nb_bias,...
    exp2a.std_obs_ratio, exp2a.obs_ratio_gm, 0, rat_changes);

%% we already plotted the lines for all refined fits for Expt 2a above, but we need the corresponding RMSEs

%the MA weight is for the obstacle obstructed target
%the PO weight is for the obstacle avoidance

%first try 66% for obstacle-obstructed target, and 66% for OA
figure; hold on; title('66% weighting');
[~, pmodel66_expt2a, ~, ~, ~, ~, ~, ~, exp2a.po.(cfld), exp2a.ma.(cfld)] = ...
    display_PO_diff_exp2a_vref( exp2a.sm_diff, exp2a.dt_diff, -1, phiA, [], 0.66, 0.66, exp2a.obs_nb_bias, exp2a.std_obs_ratio,...
    exp2a.obs_ratio_gm, {'Expt 2a: 66% weighting'});
hold off;

%try 33% for obstacle-obstructed target, and 33% for OA
figure; hold on; title('33% weighting');
[~, pmodel33_expt2a, ~, ~, ~, ~, ~, ~, exp2a.po.(cfld), exp2a.ma.(cfld)] = ...
    display_PO_diff_exp2a_vref( exp2a.sm_diff, exp2a.dt_diff, -1, phiA, [], 0.33, 0.33, exp2a.obs_nb_bias,...
    exp2a.std_obs_ratio, exp2a.obs_ratio_gm, {'Expt 2a: 33% weighting'});
hold off;

%% try cross-validation to determine the model weights

% sd1_all = [std_st_test_obs_present_comb;exp2a.s1];
% sd2_all = [std_dt_test_sub_comb;exp2a.s2];
% 
% sm_diff_all = [sm_diff; -exp2a.sm_diff];
% dt_diff_all = [dt_diff; exp2a.dt_diff];

%obs_nb_bias_all = [obs_nb_bias.comb; exp2a.obs_nb_bias.comb];

%dt_diff_ma_shifted, dt_diff_po_shifted

%% plot the Expt 2a version of Fig 4b

ma_gm.exp2a = mean(exp2a.ma_pred.sub);
ma_err.exp2a = std(exp2a.ma_pred.sub)/sqrt(exp2a.num_subjects)*1.96;

po_gm.exp2a = mean(exp2a.po_pred.sub);
po_err.exp2a = std(exp2a.po_pred.sub)/sqrt(exp2a.num_subjects)*1.96;

figure; hold on;
plot([1], dtt_obs_avg.expt2a, 's', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'markersize', 8);
display_xy_error_V2([1],dtt_obs_avg.expt2a, [], dtt_obs_err.expt2a, 'k');

%plot MA and PO as lines
hline(po_gm.exp2a, 'g');
hline(ma_gm.exp2a, 'color', [0.5,0,0.5]);

set(gca, 'XTick', [1],'XTickLabel', { 'obstacle-present 2-target trial data'});
ylabel('IMD (deg)');
title('Exp 2a data');
ylim([-15,15]);
set(gca, 'YTick', [-12:2:12],'YTickLabel', {'-12', '', '', '-6', '', '', '0', '', '', '6', '', '', '12'});
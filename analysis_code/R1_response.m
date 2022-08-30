%% plot v1 vs v2

v1.b = std_st_test_obs_present_comb.^2;
v2.b = std_dt_test_sub_comb.^2;

v1.a = exp2a.s1.^2;
v2.a = exp2a.s2.^2;

%pool the variance as well
v1.p = [v1.a; v1.b];
v2.p = [v2.a; v2.b];
total_subjects = length(v1.a) + length(v1.b);

figure; hold on;

%plot Expt 2b data
plot(sqrt(median(v1.b)), sqrt(median(v2.b)), 'bs', 'markersize', 10, 'linewidth', 1.5, 'HandleVisibility', 'Off');
plot(sqrt(mean(v1.b)), sqrt(mean(v2.b)), 'b^', 'markersize', 8, 'linewidth', 1.5, 'HandleVisibility', 'Off');
plot(sqrt(v1.b),sqrt(v2.b), 'bo', 'HandleVisibility', 'Off');

%repeat for Expt 2a
plot(sqrt(median(v1.a)), sqrt(median(v2.a)), 'rs', 'markersize', 10, 'linewidth', 1.5, 'HandleVisibility', 'Off');
plot(sqrt(mean(v1.a)), sqrt(mean(v2.a)), 'r^', 'markersize', 8, 'linewidth', 1.5, 'HandleVisibility', 'Off');
plot(sqrt(v1.a),sqrt(v2.a), 'ro', 'HandleVisibility', 'Off');

xx = [0:0.1:12];
yy = [0:0.1:12];
plot(nan,nan, 'ks','displayname', 'median');
plot(nan,nan, 'k^','displayname', 'mean');
plot(nan,nan, 'ko','displayname', 'individual subject data');
% plot(nan,nan, 'b-','displayname', 'Expt 2b');
% plot(nan,nan, 'r-','displayname', 'Expt 2a');

plot(xx,yy, 'k--','displayname', 'y=x');

leg1 = legend('show');
set(leg1, 'Location', 'Best');

% title('Expt 2b data');
title('Fig R1');

xlabel('Standard deviation during obstacle-obstructed 1-target trials');
ylabel('Standard deviation during obstacle-present 2-target trials');

%do a paired t-test for the pooled data (SD)
[~,p_var] = ttest(sqrt(v1.p), sqrt(v2.p));

%do t-test for variance
[~,p_var2] = ttest(v1.p, v2.p);

%do t-test for IQR

%calculate mean based on mean of all ratios
mean_ratio_all = mean(sqrt(v2.p)./sqrt(v1.p));

%look at p-value at an individual subject level
p_var_sub.b = nan(num_subjects,1);
p_var_sub.a = nan(exp2a.num_subjects,1);

sm_comb_all.b = (-s_AL + s_AR)/2;
dt_comb_all.b = (-left_obsA_dtt_iang_cut + right_obsA_dtt_iang_cut)/2;

for k=1:num_subjects
    [~,p_var_sub.b(k)] = vartest2(dt_comb_all.b(k,:),sm_comb_all.b(k,:), 'alpha', 0.01);
end

sm_comb_all.a = exp2a.st_all.diff;
dt_comb_all.a = exp2a.dt_all.diff;
for k=1:exp2a.num_subjects
    [~,p_var_sub.a(k)] = vartest2(dt_comb_all.a(k,:),sm_comb_all.a(k,:), 'alpha', 0.01);
end

p_var_sub.p = [p_var_sub.a; p_var_sub.b];
figure; hold on;
hist(p_var_sub.p,[0:0.01:1]);
vline(0.01, 'k', 'linewidth', 1.5);

%% Fig R1: repeat Fig 4b, but drop correction factor

%first is Expt 2a
ma_gm.exp2a = mean(exp2a.ma_pred.sub);
ma_err.exp2a = std(exp2a.ma_pred.sub)/sqrt(exp2a.num_subjects)*1.96;

po_gm.exp2a = mean(exp2a.po_pred.sub);
po_err.exp2a = std(exp2a.po_pred.sub)/sqrt(exp2a.num_subjects)*1.96;

figure; hold on;
plot([1], dtt_obs_avg.expt2a, 's', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'markersize', 6);
display_xy_error_V2([1],dtt_obs_avg.expt2a, [], dtt_obs_err.expt2a, 'k');

%plot MA and PO as lines
hline(po_gm.exp2a, 'g');
hline(ma_gm.exp2a, 'color', [0.5,0,0.5]);

%determine MA and PO without correction factor
po_ncf.exp2a = determine_po_pred_expt2a(exp2a.sm_diff, phiA, exp2a.std_obs_ratio*0+1, exp2a.obs_ratio_gm*0+1, 0.5);

%plot it as a dashed green line
hline(mean(po_ncf.exp2a.sub), 'g', 'linestyle', '--');

set(gca, 'XTick', [1],'XTickLabel', { 'obstacle-present 2-target trial data'});
ylabel('IMD (deg)');
title('Exp 2a data');
ylim([-12,20]);
set(gca, 'YTick', [-12:4:20],'YTickLabel', {'-12', '', '', '0', '', '', '12', '', '20'});

%calculate new pairwise distance
po_dist_ncf.exp2a = sqrt((exp2a.dt_diff - po_ncf.exp2a.gm).^2); %note that gm is equivalent to sub
[~, p_ncf.exp2a] = ttest(ma_dist.exp2a, po_dist_ncf.exp2a);

%% repeat above for Expt 2b

figure; hold on;
plot([1], dtt_obs_avg.expt2b, 's', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'markersize', 6);
display_xy_error_V2([1],dtt_obs_avg.expt2b, [], dtt_obs_err.expt2b, 'k');

%plot MA and PO as lines
hline(po_gm.exp2b, 'g');
hline(ma_gm.exp2b, 'color', [0.5,0,0.5]);

%determine PO without correction factor
po_ncf.exp2b = determine_po_pred(sm_diff, phiA, obs_ratio_sub*0+1, obs_ratio_gm*0+1, 0.5);

%plot it as a dashed green line
hline(mean(po_ncf.exp2b.sub), 'g', 'linestyle', '--');

set(gca, 'XTick', [1],'XTickLabel', { 'obstacle-present 2-target trial data'});
ylabel('IMD (deg)');
title('Exp 2b data');
ylim([-12,20]);
set(gca, 'YTick', [-12:4:20],'YTickLabel', {'-12', '', '', '0', '', '', '12', '', '20'});

%calculate new pairwise distance
po_dist_ncf.exp2b = sqrt((dt_diff - po_ncf.exp2b.gm).^2);
[~, p_ncf.exp2b] = ttest(ma_dist.exp2b, po_dist_ncf.exp2b);

%% make bar plots of mean and SEM of of motor variability for pooled, 2b, and 2a data for 1 vs 2-target trials
%get corresponding p values as well

%calculate the mean and SEM of the pooled data
mean_v1p = mean(sqrt(v1.p));
sem_v1p = std(sqrt(v1.p))/sqrt(total_subjects);

mean_v2p = mean(sqrt(v2.p));
sem_v2p = std(sqrt(v2.p))/sqrt(total_subjects);

%repeat for Expt 2a
mean_v1a = mean(sqrt(v1.a));
sem_v1a = std(sqrt(v1.a))/sqrt(exp2a.num_subjects);

mean_v2a = mean(sqrt(v2.a));
sem_v2a = std(sqrt(v2.a))/sqrt(exp2a.num_subjects);

%repeat for Expt 2b
mean_v1b = mean(sqrt(v1.b));
sem_v1b = std(sqrt(v1.b))/sqrt(num_subjects);

mean_v2b = mean(sqrt(v2.b));
sem_v2b = std(sqrt(v2.b))/sqrt(num_subjects);

figure; hold on;
title('Pooled data');

hb1 = barwitherr( [sem_v1p, sem_v2p],...
    [mean_v1p, mean_v2p ] );

set(gca, 'XTick', [1,2],'XTickLabel', {'obstacle-obstructed 1TT', 'obstacle-present 2TT'});
ylabel('SD (deg)');
ylim([0,8]);

%get p value
[~,p_vp] = ttest(sqrt(v1.p), sqrt(v2.p));

%%%try Expt 2a
figure; hold on;
title('Expt 2a data');

hb1 = barwitherr( [sem_v1a, sem_v2a],...
    [mean_v1a, mean_v2a ] );

set(gca, 'XTick', [1,2],'XTickLabel', {'obstacle-obstructed 1TT', 'obstacle-present 2TT'});
ylabel('SD (deg)');
ylim([0,8]);

%get p value
[~,p_va] = ttest(sqrt(v1.a), sqrt(v2.a));

%%%try Expt 2b
figure; hold on;
title('Expt 2b data');

hb1 = barwitherr( [sem_v1b, sem_v2b],...
    [mean_v1b, mean_v2b ] );

set(gca, 'XTick', [1,2],'XTickLabel', {'obstacle-obstructed 1TT', 'obstacle-present 2TT'});
ylabel('SD (deg)');
ylim([0,8]);

%get p value
[~,p_vb] = ttest(sqrt(v1.b), sqrt(v2.b));

%% for each subject, generate a histogram of p-values corresponding to the difference between left and right
% obstacle-obstructed 1-target trials

s_stL.b = sqrt(var_st_test_obs_present(1,:));
s_stR.b = sqrt(var_st_test_obs_present(2,:));

%calculate left/right variance for Expt 2a
s_stL.a = nanstd(exp2a.st_all.L,0,2);
s_stR.a = nanstd(exp2a.st_all.R,0,2);

s_stL.p = [s_stL.a(:);s_stL.b(:)];
s_stR.p = [s_stR.a(:);s_stR.b(:)];

%calculate percent difference for each person
lvr_diff = abs(s_stL.p-s_stR.p)./ (s_stL.p/2+s_stR.p/2) * 100;

p_var_st_sub.b = nan(num_subjects,1);
p_var_st_sub.a = nan(exp2a.num_subjects,1);

for k=1:num_subjects
    [~,p_var_st_sub.b(k)] = vartest2(-s_AL(k,:),s_AR(k,:), 'alpha', 0.01);
end

for k=1:exp2a.num_subjects
    [~,p_var_st_sub.a(k)] = vartest2(-exp2a.st_all.L(k,:),exp2a.st_all.R(k,:), 'alpha', 0.01);
end

p_var_st_sub.p = [p_var_st_sub.a; p_var_st_sub.b];
figure; hold on;
hist(p_var_st_sub.p,[0:0.01:1]);
vline(0.01, 'k', 'linewidth', 1.5);

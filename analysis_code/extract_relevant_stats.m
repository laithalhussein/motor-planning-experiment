close all;
home;

%data is subject x sample

%% get mean and se for all conditions

%left obstacle data
left_tgt_left_obs = [mean(nanmean(iang_st_left_obsA_bl,2)), mean(nanstd(iang_st_left_obsA_bl,0,1))/sqrt(num_subjects)];
right_tgt_left_obs = [mean(nanmean(iang_st_right_obsA_nb,2)), mean(nanstd(iang_st_right_obsA_nb,0,1))/sqrt(num_subjects)];
cntr_tgt_left_obs = [mean(nanmean(iang_strt_obsA_left,2)), mean(nanstd(iang_strt_obsA_left,0,1))/sqrt(num_subjects)];

%right obstacle data
left_tgt_right_obs = [mean(nanmean(iang_st_left_obsA_nb,2)), mean(nanstd(iang_st_left_obsA_nb,0,1))/sqrt(num_subjects)];
right_tgt_right_obs = [mean(nanmean(iang_st_right_obsA_bl,2)), mean(nanstd(iang_st_right_obsA_bl,0,1))/sqrt(num_subjects)];
cntr_tgt_right_obs = [mean(nanmean(iang_strt_obsA_right,2)), mean(nanstd(iang_strt_obsA_right,0,1))/sqrt(num_subjects)];

%DTT
dtt_left_obs = [mean(nanmean([iang_dt_left_obsA_bl, iang_dt_right_obsA_nb],2)), mean(nanstd([iang_dt_left_obsA_bl, iang_dt_right_obsA_nb],0,1))/sqrt(num_subjects)];
dtt_right_obs = [mean(nanmean([iang_dt_right_obsA_bl, iang_dt_left_obsA_nb],2)), mean(nanstd([iang_dt_right_obsA_bl, iang_dt_left_obsA_nb],0,1))/sqrt(num_subjects)];

%baseline data
left_tgt_bln = [mean(nanmean(iang_st_left_bln,2)), mean(nanstd(iang_st_left_bln,0,1))/sqrt(num_subjects)];
right_tgt_bln = [mean(nanmean(iang_st_right_bln,2)), mean(nanstd(iang_st_right_bln,0,1))/sqrt(num_subjects)];
cntr_tgt_bln = [mean(nanmean(iang_st_strt_bln,2)), mean(nanstd(iang_st_strt_bln,0,1))/sqrt(num_subjects)];

dtt_bln = [mean(nanmean(iang_dt_bln,2)), mean(nanstd(iang_dt_bln,0,1))/sqrt(num_subjects)];


%% perform t-tests

%t-test between all 2 tgt trials and 0
[h1, p1] = ttest([nanmean([iang_dt_left_obsA_bl, iang_dt_right_obsA_nb, iang_dt_right_obsA_nb, iang_dt_left_obsA_nb],2)]);

%repeat to test left obs against baseline
[h2, p2] = ttest2([nanmean([iang_dt_left_obsA_bl, iang_dt_right_obsA_nb],2)], nanmean(iang_dt_bln,2),'Vartype','unequal');

%repeat to test rigt obs against baseline
[h11, p11] = ttest2([nanmean([iang_dt_right_obsA_bl, iang_dt_left_obsA_nb],2)], nanmean(iang_dt_bln,2),'Vartype','unequal');

%repeat and test bln 2tt against bln cntr data
[h10, p10] = ttest2(nanmean(iang_st_strt_bln,2), nanmean(iang_dt_bln,2),'Vartype','unequal');

%test each tgt's data from obstacle to baseline, and do it separately for left vs right obstacle

%left obstacle, left tgt
[h3, p3] = ttest2(nanmean(iang_st_left_obsA_bl,2), nanmean(iang_st_left_bln,2),'Vartype','unequal');

%left obstacle, right tgt
[h4, p4] = ttest2(nanmean(iang_st_right_obsA_nb,2), nanmean(iang_st_right_bln,2));

%left obstacle, cntr tgt
[h5, p5] = ttest2(nanmean(iang_strt_obsA_left,2), nanmean(iang_st_strt_bln,2));

%right obstacle, left tgt
[h6, p6] = ttest2(nanmean(iang_st_left_obsA_nb,2), nanmean(iang_st_left_bln,2));

%right obstacle, right tgt
[h7, p7] = ttest2(nanmean(iang_st_right_obsA_bl,2), nanmean(iang_st_right_bln,2),'Vartype','unequal');

%right obstacle, cntr tgt
[h8, p8] = ttest2(nanmean(iang_strt_obsA_right,2), nanmean(iang_st_strt_bln,2));

%% permutation test between strt 1-tgt trials (bln) and bln 2-target trials

p_var_trial_type = nan(num_subjects,1);
pse_var_trial_type = nan(num_subjects, 1);

for kq2 = 1:num_subjects
    [p_var_trial_type(kq2), pse_var_trial_type(kq2)] = permtest_var_1210_2019a(iang_st_strt_bln(kq2,:), iang_dt_bln(kq2,:), 1);
    %[~,p_var_trial_type(kq2)] = vartest2(right_obsA_dtt_iang(kq2,:), left_obsA_dtt_iang(kq2,:));  
end

%% try pooling the data before testing

a = iang_st_strt_bln(:);
b = iang_dt_bln(:);

[p9, p9se] = permtest_var_1210_2019a(a, b, 0, 1e4);

%%
% keyboard;



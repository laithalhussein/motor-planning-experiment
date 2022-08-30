close all;
home;

%test specific distributions to see if they're normal
%use One-sample Kolmogorov-Smirnov test, input must be standardized

%only need to test data that are used for statistical comparisons or are
%thrown into a linear regression

% [h1,p1] = kstest(sm_diff-nanmean(sm_diff)./nanstd(sm_diff));

%obstacle-present 2-target trials
[h1,p1] = kstest(dt_diff-nanmean(dt_diff)./nanstd(dt_diff)); %fails

%obstacle-free 2-target trials
of2 = nanmedian(iang_dt_bln,2);
[h2,p2] = kstest((of2-mean(of2))./std(of2)); %passes

%RMSEs
[h3,p3] = kstest((ma_dist.exp2b-mean(ma_dist.exp2b))./std(ma_dist.exp2b)); %passes
[h4,p4] = kstest((po_dist.exp2b-mean(po_dist.exp2b))./std(po_dist.exp2b)); %passes

%% repeat for Expt 2a

%obstacle-present 2-target trials
[h5,p5] = kstest(exp2a.dt_diff-nanmean(exp2a.dt_diff)./nanstd(exp2a.dt_diff)); %fails

%obstacle-free 2-target trials
of2_exp2a = nanmedian(exp2a.iang_dt_bln,2);
[h6,p6] = kstest((of2_exp2a-mean(of2_exp2a))./std(of2_exp2a)); %passes

%RMSEs
[h7,p7] = kstest((ma_dist.exp2a-mean(ma_dist.exp2a))./std(ma_dist.exp2a)); %passes
[h8,p8] = kstest((po_dist.exp2a-mean(po_dist.exp2a))./std(po_dist.exp2a)); %passes

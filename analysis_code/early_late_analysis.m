%% early vs late PO and MA comparison 
%NOTE: we might want to repeat a lot of the anlyses in initial_movement_analysis for a specific period

% do the predictions we have above, but this time only use early and late samples

early_frac = 0.25; %fraction of trials to use from the beginning
late_frac = 0.25; %same for the late period

%we could divide the data into quarters to see if that has an effect
%idx_st_all = reshape([1:size(iang_st_left_obsA_bl,2)], size(iang_st_left_obsA_bl,2)/4, 4);

%early_period_obsA_idx = idx_st_all(:,3);
%late_period_obsA_idx = idx_st_all(:,4);

 early_period_obsA_idx = [1: floor(length(iang_st_left_obsA_bl) * early_frac)];
 %early_period_obsB_idx = [1: floor(length(iang_st_left_obsB_bl) * early_frac)]; %last quarter for these are only 2 trials!

late_period_obsA_idx = [length(iang_st_left_obsA_bl) - floor(length(iang_st_left_obsA_bl) * late_frac) : length(iang_st_left_obsA_bl)];
%late_period_obsB_idx = [length(iang_st_left_obsB_bl) - floor(length(iang_st_left_obsB_bl) * late_frac) : length(iang_st_left_obsB_bl)];

early_dat.st_left_A = st_block_dat.obsA_L(:,early_period_obsA_idx);
early_dat.st_right_A = st_block_dat.obsA_R(:,early_period_obsA_idx);
% early_dat.st_left_B = st_block_dat.obsB_L(:,early_period_obsB_idx);
% early_dat.st_right_B = st_block_dat.obsB_R(:,early_period_obsB_idx);

early_dat.dt_left_A = left_obsA_dtt_iang(:, early_period_obsA_idx); %since there are more dt trials, this is not actually using the last quarter!
early_dat.dt_right_A = right_obsA_dtt_iang(:, early_period_obsA_idx);
% early_dat.dt_left_B = left_obsB_dtt_iang(:, early_period_obsB_idx);
% early_dat.dt_right_B = right_obsB_dtt_iang(:, early_period_obsB_idx);

late_dat.st_left_A = st_block_dat.obsA_L(:,late_period_obsA_idx);
late_dat.st_right_A = st_block_dat.obsA_R(:,late_period_obsA_idx);
% late_dat.st_left_B = st_block_dat.obsB_L(:,late_period_obsB_idx);
% late_dat.st_right_B = st_block_dat.obsB_R(:,late_period_obsB_idx);

late_dat.dt_left_A = left_obsA_dtt_iang(:, late_period_obsA_idx);
late_dat.dt_right_A = right_obsA_dtt_iang(:, late_period_obsA_idx);
% late_dat.dt_left_B = left_obsB_dtt_iang(:, late_period_obsB_idx);
% late_dat.dt_right_B = right_obsB_dtt_iang(:, late_period_obsB_idx);


%%%%%%%%%%% plot them
figure; hold on; title('Early vs Late: Obs A');

%do it for obs A first
display_optimal_prediction_EL( early_dat.st_left_A, early_dat.dt_left_A, late_dat.st_left_A, late_dat.dt_left_A, -1, phiA, [] );
display_optimal_prediction_EL( early_dat.st_right_A, early_dat.dt_right_A, late_dat.st_right_A, late_dat.dt_right_A, 1, phiA, [] );
% 
% %now for obstacle B
% figure; hold on; title('Early vs Late: Obs B');
% 
% display_optimal_prediction_EL( early_dat.st_left_B, early_dat.dt_left_B, late_dat.st_left_B, late_dat.dt_left_B, -1, phiB, [] );
% display_optimal_prediction_EL( early_dat.st_right_B, early_dat.dt_right_B, late_dat.st_right_B, late_dat.dt_right_B, 1, phiB, [] );
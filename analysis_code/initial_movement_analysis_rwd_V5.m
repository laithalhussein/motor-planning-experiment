%additional analysis for the obstacle stuff
%run the other analysis before this one!!

close all;

clearvars -except Cdr dat_all exp_seq Ideal_dir_all idx_all info_all MT_labels_sub_all...
    num_subjects optimal_indicator_all Pos_sub_all S Sub_id_all xTraj_sub_all yTraj_sub_all...
    rwd_label

red_=[1 0 0];
purple=[128,0,128]/256;
violet = [128, 130, 238] / 256;
grey=[0.5 0.5 0.5];
orange=[255,140,0]/256;
black_=[0,0,0];
dblue = [0,0,255] / 256;
cyan = [0,255,255] / 256;

%calculate the angle for all the cases at a certain point into the movement

initial_marker = 4; %cm into the movement

%bias_flag = 0; %do we want to subtract by individual subject biases?
st_bias_flag = 1;
dt_bias_flag = 0; 

dt_bias_type = 2; %set this to 2 if we want to use the null trials during the test period

inside_filter = 1; %option to remove trials where subjects move inside relative to the obstacle on STT blocked trials

%assuming that all subjects completed the same number of trials
num_trials = length(Pos_sub_all(1,:));

%calculate the angle and save it
initial_angle = nan(num_subjects, num_trials);
%v1 = [0,initial_marker]; %vector we will calculate the angle wrt

initial_mov = cell(num_subjects, num_trials);

for qq = 1:num_subjects
    for kn = 1:num_trials
        
        current_pos = Pos_sub_all{qq,kn};
        initial_idx = find(current_pos >= initial_marker,1,'first');
        %now use this to get the vector for angle calculations
        x_data = (xTraj_sub_all{qq,kn});
        y_data = yTraj_sub_all{qq,kn};
        
        v1 = [x_data(1) - S(1),initial_marker]; %vector we will calculate the angle wrt
        v2 = [x_data(initial_idx) - S(1), y_data(initial_idx) - S(2)]; %vector corresponding to movement
        
%         if kn == 10, keyboard; end
        
        if ~isempty(v2)
            angle = atan2d(det([v1;v2]),dot(v1,v2)); %CCCW is positive, so multiply by -1 (??)
        else
%             keyboard; %these are bad movements
            angle= nan;
        end
        
        initial_angle(qq,kn) = angle;
       
        initial_mov{qq,kn}(:,1) = x_data(1:initial_idx);
        initial_mov{qq,kn}(:,2) = y_data(1:initial_idx);
                
    end
end

%% since its used a lot, make a matrix that contains the rewarded trials for each subject

sub_rwd_all = nan(num_subjects, size(MT_labels_sub_all,2) )';

for qn = 1:num_subjects
    
    tmp_bad = squeeze(info_all.Bad_labels_sub_all(qn,:,:));
    
    %lets not count glass / horn for now
    %tmp_bad = tmp_bad(:,1) | tmp_bad(:,2) | tmp_bad(:,3); %combine them
    tmp_bad = tmp_bad(:,3); %only looking at obstacle hitting
    
    
    MT_label_tmp = MT_labels_sub_all(qn,:)'; % 1 means it was too slow
    MT = dat_all.MT_sub_all(qn,:);
    
    TF = find(MT <= .225); %too fast
    MT_label_tmp(TF) = 1; %also bad
    
    sub_rwd_all(:,qn) = MT_label_tmp | tmp_bad;
    
end

%% also need to average baseline trajectories!!!
%get the indices

%NOTE: we could use the above to find the average angle between trials, but
%this is more suceptible to noise, so lets average the baseline
%trajectories for each subject, and then calculate that angle

%IMPORTANT: probably need to interpolate

xbln_st_left = xTraj_sub_all(:, idx_all.st_left_bln );
ybln_st_left = yTraj_sub_all(:, idx_all.st_left_bln );

xbln_st_right = xTraj_sub_all(:, idx_all.st_right_bln );
ybln_st_right = yTraj_sub_all(:, idx_all.st_right_bln );

xbln_st_strt = xTraj_sub_all(:, idx_all.st_strt_bln );
ybln_st_strt = yTraj_sub_all(:, idx_all.st_strt_bln );

xbln_dt = xTraj_sub_all(:, idx_all.dt_bln );
ybln_dt = yTraj_sub_all(:, idx_all.dt_bln );

xtest_null_dt = xTraj_sub_all(:, idx_all.dt_null_test );
ytest_null_dt = yTraj_sub_all(:, idx_all.dt_null_test );

%we will average b/w trials for each subject to calculate the initial
%angle, so lets take only a certain amount of samples b/w trials

max_samples = 95; %this should be more than enough for the inital angle

xbln_st_left_tmp = cellfun(@(x) x(1:max_samples), xbln_st_left,'UniformOutput',false);
ybln_st_left_tmp = cellfun(@(x) x(1:max_samples), ybln_st_left,'UniformOutput',false);
xbln_st_right_tmp = cellfun(@(x) x(1:max_samples), xbln_st_right,'UniformOutput',false);
ybln_st_right_tmp = cellfun(@(x) x(1:max_samples), ybln_st_right,'UniformOutput',false);
xbln_st_strt_tmp = cellfun(@(x) x(1:max_samples), xbln_st_strt,'UniformOutput',false);
ybln_st_strt_tmp = cellfun(@(x) x(1:max_samples), ybln_st_strt,'UniformOutput',false);
xbln_dt_tmp = cellfun(@(x) x(1:max_samples), xbln_dt,'UniformOutput',false);
ybln_dt_tmp = cellfun(@(x) x(1:max_samples), ybln_dt,'UniformOutput',false);

xtest_null_dt_tmp = cellfun(@(x) x(1:max_samples), xtest_null_dt,'UniformOutput',false);
ytest_null_dt_tmp = cellfun(@(x) x(1:max_samples), ytest_null_dt,'UniformOutput',false);

%preallocation
xbln_st_left_avg = nan(num_subjects, max_samples);
ybln_st_left_avg = nan(num_subjects, max_samples);
xbln_st_right_avg = nan(num_subjects, max_samples);
ybln_st_right_avg = nan(num_subjects, max_samples);

xbln_st_strt_avg = nan(num_subjects, max_samples);
ybln_st_strt_avg = nan(num_subjects, max_samples);

xbln_dt_avg = nan(num_subjects, max_samples);
ybln_dt_avg = nan(num_subjects, max_samples);

xtest_dt_null_avg = nan(num_subjects, max_samples);
ytest_dt_null_avg = nan(num_subjects, max_samples);

iang_st_left_bln = nan(num_subjects, 1);
iang_st_right_bln = nan(num_subjects, 1);
iang_st_strt_bln = nan(num_subjects, 1);
iang_dt_bln = nan(num_subjects, 1);
iang_dt_null_test = nan(num_subjects, 1);

for jk = 1:num_subjects
%    xbln_st_left_avg(jk,:) =  nanmean( cell2mat(xbln_st_left_tmp(jk,:)), 2);
%    ybln_st_left_avg(jk,:) =  nanmean( cell2mat(ybln_st_left_tmp(jk,:)), 2);
%    
%    xbln_st_right_avg(jk,:) =  nanmean( cell2mat(xbln_st_right_tmp(jk,:)), 2);
%    ybln_st_right_avg(jk,:) =  nanmean( cell2mat(ybln_st_right_tmp(jk,:)), 2);
%    
%    xbln_st_strt_avg(jk,:) =  nanmean( cell2mat(xbln_st_strt_tmp(jk,:)), 2);
%    ybln_st_strt_avg(jk,:) =  nanmean( cell2mat(ybln_st_strt_tmp(jk,:)), 2);
%    
%    xbln_dt_avg(jk,:) =  nanmean( cell2mat(xbln_dt_tmp(jk,:)), 2);
%    ybln_dt_avg(jk,:) =  nanmean( cell2mat(ybln_dt_tmp(jk,:)), 2);


    xbln_st_left_avg(jk,:) =  nanmedian( cell2mat(xbln_st_left_tmp(jk,:)), 2);
    ybln_st_left_avg(jk,:) =  nanmedian( cell2mat(ybln_st_left_tmp(jk,:)), 2);

    xbln_st_right_avg(jk,:) =  nanmedian( cell2mat(xbln_st_right_tmp(jk,:)), 2);
    ybln_st_right_avg(jk,:) =  nanmedian( cell2mat(ybln_st_right_tmp(jk,:)), 2);

    xbln_st_strt_avg(jk,:) =  nanmedian( cell2mat(xbln_st_strt_tmp(jk,:)), 2);
    ybln_st_strt_avg(jk,:) =  nanmedian( cell2mat(ybln_st_strt_tmp(jk,:)), 2);

    xbln_dt_avg(jk,:) =  nanmedian( cell2mat(xbln_dt_tmp(jk,:)), 2);
    ybln_dt_avg(jk,:) =  nanmedian( cell2mat(ybln_dt_tmp(jk,:)), 2);
    
    xtest_dt_null_avg(jk,:) =  nanmedian( cell2mat(xtest_null_dt_tmp(jk,:)), 2);
    ytest_dt_null_avg(jk,:) =  nanmedian( cell2mat(ytest_null_dt_tmp(jk,:)), 2);
   
   %calculate the angle based on the average trajectory
   P1 = cumsum(sqrt(diff(xbln_st_left_avg(jk,:)).^2+diff(ybln_st_left_avg(jk,:)).^2));
   P2 = cumsum(sqrt(diff(xbln_st_right_avg(jk,:)).^2+diff(ybln_st_right_avg(jk,:)).^2));
   P3 = cumsum(sqrt(diff(xbln_st_strt_avg(jk,:)).^2+diff(ybln_st_strt_avg(jk,:)).^2));
   P4 = cumsum(sqrt(diff(xbln_dt_avg(jk,:)).^2+diff(ybln_dt_avg(jk,:)).^2));
   P5 = cumsum(sqrt(diff(xtest_dt_null_avg(jk,:)).^2+diff(ytest_dt_null_avg(jk,:)).^2));
   
   idx1 = find(P1 >=initial_marker,1,'first');
   idx2 = find(P2 >=initial_marker,1,'first');
   idx3 = find(P3 >=initial_marker,1,'first');
   idx4 = find(P4 >=initial_marker,1,'first');
   idx5 = find(P5 >=initial_marker,1,'first');
   
   v2_st_left = [ xbln_st_left_avg(jk,idx1) - S(1), ybln_st_left_avg(jk,idx1) - S(2) ];
   v2_st_right = [ xbln_st_right_avg(jk,idx2) - S(1), ybln_st_right_avg(jk,idx2) - S(2) ];
   v2_st_strt = [ xbln_st_strt_avg(jk,idx3) - S(1), ybln_st_strt_avg(jk,idx3) - S(2) ];
   v2_dt = [ xbln_dt_avg(jk,idx4) - S(1), ybln_dt_avg(jk,idx4) - S(2) ];
   v2_dt_null_test = [ xtest_dt_null_avg(jk,idx5) - S(1), ytest_dt_null_avg(jk,idx5) - S(2) ];
   
   %get the angle
   iang_st_left_bln(jk) = atan2d(det([v1;v2_st_left]),dot(v1,v2_st_left));
   iang_st_right_bln(jk) = atan2d(det([v1;v2_st_right]),dot(v1,v2_st_right));
   iang_st_strt_bln(jk) = atan2d(det([v1;v2_st_strt]),dot(v1,v2_st_strt));
   iang_dt_bln(jk) = atan2d(det([v1;v2_dt]),dot(v1,v2_dt));
   iang_dt_null_test(jk) = atan2d(det([v1;v2_dt_null_test]),dot(v1,v2_dt_null_test));
   
end
        
%now we can calculate the bias...
%using the average trajectory might not be as robust as I think (without interpolating!)
bias_st_left = (iang_st_left_bln - 30) * 1;
bias_st_right = (iang_st_right_bln + 30) * 1;
bias_st_strt = iang_st_strt_bln + 0;
bias_dt_null_bln = iang_dt_bln + 0;
bias_dt_null_test = iang_dt_null_test + 0;

if dt_bias_type ==1
    bias_dt_null = bias_dt_null_bln;
elseif dt_bias_type==2
    bias_dt_null = bias_dt_null_test;
end


%% now index these for all our cases

% first for the angles
iang_st_left_obsA_bl = initial_angle(:,idx_all.st_left_obsA_bl);
iang_st_right_obsA_bl = initial_angle(:,idx_all.st_right_obsA_bl);
iang_st_left_obsB_bl = initial_angle(:,idx_all.st_left_obsB_bl);
iang_st_right_obsB_bl = initial_angle(:,idx_all.st_right_obsB_bl);

iang_st_left_obsA_nb = initial_angle(:,idx_all.st_left_obsA_nb);
iang_st_right_obsA_nb = initial_angle(:,idx_all.st_right_obsA_nb);
iang_st_left_obsB_nb = initial_angle(:,idx_all.st_left_obsB_nb);
iang_st_right_obsB_nb = initial_angle(:,idx_all.st_right_obsB_nb);

iang_dt_left_obsA_bl = initial_angle(:,idx_all.dt_left_obsA_bl);
iang_dt_right_obsA_bl = initial_angle(:,idx_all.dt_right_obsA_bl);
iang_dt_left_obsB_bl = initial_angle(:,idx_all.dt_left_obsB_bl);
iang_dt_right_obsB_bl = initial_angle(:,idx_all.dt_right_obsB_bl);

iang_dt_left_obsA_nb = initial_angle(:,idx_all.dt_left_obsA_nb);
iang_dt_right_obsA_nb = initial_angle(:,idx_all.dt_right_obsA_nb);
iang_dt_left_obsB_nb = initial_angle(:,idx_all.dt_left_obsB_nb);
iang_dt_right_obsB_nb = initial_angle(:,idx_all.dt_right_obsB_nb);

iang_strt_obsA_left = initial_angle(:,idx_all.strt_obsA_left);
iang_strt_obsA_right = initial_angle(:,idx_all.strt_obsA_right);
iang_strt_obsB_left = initial_angle(:,idx_all.strt_obsB_left);
iang_strt_obsB_right = initial_angle(:,idx_all.strt_obsB_right);

fn = 210; %figure number for this analysis to start from
figure(fn);
figure(fn+1);
figure(fn+2);
figure(fn+3);

c1 = 0; c2 = 1; c3 = 1;

MT_opt_obsA = [];       MT_opt_obsB = [];
MT_ma_obsA = [];       MT_ma_obsB = [];

p_dtt_null_mt_all = nan(num_subjects,1);
p_dtt_A_mt_all = nan(num_subjects,1);
p_dtt_B_mt_all = nan(num_subjects,1);

xb1 = {['All'], ['STT obs A (blocked)'], ['STT obs B (blocked)'], ['DTT (null)'], ['DTT obs A (all)'], ['DTT obs B (all)']};

ss=1;
ss2 = 1 + num_subjects;

inside_indc.obsA_L = [];
inside_indc.obsA_R = [];
inside_indc.obsB_L = [];
inside_indc.obsB_R = [];

for kk = 1:num_subjects
    
    opt_idx = optimal_indicator_all(kk,:);
    %ma_idx = optimal_indicator_all(kk,:);
    
    opt_indc_st_left_obsA = opt_idx(idx_all.st_left_obsA_bl) ;
    opt_indc_st_right_obsA = opt_idx(idx_all.st_right_obsA_bl);
    opt_indc_st_left_obsB = opt_idx(idx_all.st_left_obsB_bl);
    opt_indc_st_right_obsB = opt_idx(idx_all.st_right_obsB_bl);
    
    opt_st_left_obsA = initial_angle(kk,idx_all.st_left_obsA_bl(opt_indc_st_left_obsA==1)) - (bias_st_left(kk) * st_bias_flag);
    opt_st_right_obsA = initial_angle(kk,idx_all.st_right_obsA_bl(opt_indc_st_right_obsA==1)) - (bias_st_right(kk) * st_bias_flag);
    opt_st_left_obsB = initial_angle(kk,idx_all.st_left_obsB_bl(opt_indc_st_left_obsB==1)) - (bias_st_left(kk) * st_bias_flag);
    opt_st_right_obsB = initial_angle(kk,idx_all.st_right_obsB_bl(opt_indc_st_right_obsB==1)) - (bias_st_right(kk) * st_bias_flag);
    
    ma_st_left_obsA = initial_angle(kk,idx_all.st_left_obsA_bl(opt_indc_st_left_obsA==0)) - (bias_st_left(kk) * st_bias_flag);
    ma_st_right_obsA = initial_angle(kk,idx_all.st_right_obsA_bl(opt_indc_st_right_obsA==0)) - (bias_st_right(kk) * st_bias_flag);
    ma_st_left_obsB = initial_angle(kk,idx_all.st_left_obsB_bl(opt_indc_st_left_obsB==0)) - (bias_st_left(kk) * st_bias_flag);
    ma_st_right_obsB = initial_angle(kk,idx_all.st_right_obsB_bl(opt_indc_st_right_obsB==0)) - (bias_st_right(kk) * st_bias_flag);
    
    %save the inside vs outside indices for STT trials, when people when
    %outside, we might not want to count this in our analysis
    inside_indc.obsA_L(:,kk) = ~opt_indc_st_left_obsA; % 1 means they didnt do the optimal (intended) movement
    inside_indc.obsA_R(:,kk) = ~opt_indc_st_right_obsA;
    inside_indc.obsB_L(:,kk) = ~opt_indc_st_left_obsB;
    inside_indc.obsB_R(:,kk) = ~opt_indc_st_right_obsB;
    
    
%     if kk==5, keyboard; end
    
    figure(fn);
    subplot(num_subjects,2,kk+c1); hold on;
    create_hist_init(opt_st_left_obsA ,violet, 30);
    create_hist_init(ma_st_left_obsA ,purple, 30);
    if (kk+c1) == 1, title('Left target'); end
    
    subplot(num_subjects,2,kk+c2); hold on;
    create_hist_init(opt_st_right_obsA ,orange, 30);
    create_hist_init(ma_st_right_obsA ,red_, 30);
    if (kk+c2) == 2, title('Right target'); end
    
    %repeat for obstacle B
    figure(fn+1);
    subplot(num_subjects,2,kk+c1); hold on;
    create_hist_init(opt_st_left_obsB ,violet, 30);
    create_hist_init(ma_st_left_obsB ,purple, 30);
    if (kk+c1) == 1, title('Left target'); end
    
    
    subplot(num_subjects,2,kk+c2); hold on;
    create_hist_init(opt_st_right_obsB ,orange, 30);
    create_hist_init(ma_st_right_obsB ,red_, 30);
    if (kk+c2) == 2, title('Right target'); end
    
    
    %now for DTTs...
    %For these cases, we will pool the data based on where the obstacle, not which tgt is cued, because presumably they dont know...
    figure(fn+2);
    subplot(num_subjects,2,kk+c1); hold on;
    left_obsA_dtt_iang(kk,:) = [iang_dt_left_obsA_bl(kk,:),iang_dt_right_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
    create_hist_init([iang_dt_left_obsA_bl(kk,:),iang_dt_right_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag) ,purple, 30);
    if (kk+c1) == 1, title('Leftward obstacle'); end
    
    subplot(num_subjects,2,kk+c2); hold on;
    right_obsA_dtt_iang(kk,:) = [iang_dt_right_obsA_bl(kk,:),iang_dt_left_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
    create_hist_init([iang_dt_right_obsA_bl(kk,:),iang_dt_left_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag),red_, 30);
    if (kk+c2) == 2, title('Rightward obstacle'); end
    
    %DTT for obstacle B
    figure(fn+3);
    subplot(num_subjects,2,kk+c1); hold on;
    left_obsB_dtt_iang(kk,:) = [iang_dt_left_obsB_bl(kk,:),iang_dt_right_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
    create_hist_init([iang_dt_left_obsB_bl(kk,:),iang_dt_right_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag),purple, 30);
    if (kk+c1) == 1, title('Leftward obstacle'); end
    
    subplot(num_subjects,2,kk+c2); hold on;
    right_obsB_dtt_iang(kk,:) = [iang_dt_right_obsB_bl(kk,:),iang_dt_left_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
    create_hist_init([iang_dt_right_obsB_bl(kk,:),iang_dt_left_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag),red_, 30);
    if (kk+c2) == 2, title('Rightward obstacle'); end
    
    
%     %make a 2 x2 scatter plot of the average of the median angles for STT
%     %left and right targets for each subject (motor averaging prediction)
%     %vs the observed angle on DTTs
%     figure(fn + 10);
%     %first is obs A Right target
%     subplot(2,2,1); hold on;
%     tmp_ang = median
    
    
    %if kk==2, keyboard; end
    
    c1 = c1 +1;
    c2 = c2+1;
    
    %now take care of MTs
    sub_bad_trial = squeeze(info_all.Bad_labels_sub_all(kk,:,:));
    sub_good_trial = ~sub_bad_trial;
    MT_label_sub = MT_labels_sub_all(kk,:)';
    MT = dat_all.MT_sub_all(kk,:);
    
    %focus only only sucessfully fast movements?
    %good_label= prod(sub_good_trial,2) .*MT_label_sub;
    good_label= prod(sub_good_trial,2);
    
    MT_good = MT'.*good_label;
    %perc_analysis = sum(MT_good) / length(MT)
    MT_good(MT_good==0) = nan;
    
    %now get the data
    MT_test = MT_good(exp_seq{3}(end)+1:end);
    MT_st_obsA_blk = MT_good( [idx_all.st_left_obsA_bl, idx_all.st_right_obsA_bl] );
    MT_st_obsB_blk = MT_good( [idx_all.st_left_obsB_bl, idx_all.st_right_obsB_bl] );
    MT_dt_null_test = MT_good( idx_all.dt_null_test);
    MT_dt_obsA = MT_good( [idx_all.dt_left_obsA_bl, idx_all.dt_right_obsA_bl, idx_all.dt_left_obsA_nb, idx_all.dt_right_obsA_nb]);
    MT_dt_obsB = MT_good( [idx_all.dt_left_obsB_bl, idx_all.dt_right_obsB_bl, idx_all.dt_left_obsB_nb, idx_all.dt_right_obsB_nb]);
    
    %also separate them based on optimal vs averaging as well, pooled across all subjects
    MT_opt_obsA_tmp = MT_st_obsA_blk( [(opt_indc_st_left_obsA==1), (opt_indc_st_right_obsA==1)]);
    MT_opt_obsB_tmp = MT_st_obsB_blk( [(opt_indc_st_left_obsB==1), (opt_indc_st_right_obsB==1)]);
    
    MT_ma_obsA_tmp = MT_st_obsA_blk( [(opt_indc_st_left_obsA==0), (opt_indc_st_right_obsA==0)]);
    MT_ma_obsB_tmp = MT_st_obsB_blk( [(opt_indc_st_left_obsB==0), (opt_indc_st_right_obsB==0)]);
    
    MT_opt_obsA = [MT_opt_obsA; MT_opt_obsA_tmp];
    MT_opt_obsB = [MT_opt_obsB; MT_opt_obsB_tmp];
    MT_ma_obsA = [MT_ma_obsA; MT_ma_obsA_tmp];
    MT_ma_obsB = [MT_ma_obsB; MT_ma_obsB_tmp];
    
    
    %separate optimal vs averaging mt based on target direction
    MT_st_right_obsA_blk = MT_good( [idx_all.st_right_obsA_bl] );
    MT_st_left_obsA_blk = MT_good( [idx_all.st_left_obsA_bl] );
    
    MT_st_right_obsB_blk = MT_good( [idx_all.st_right_obsB_bl] );
    MT_st_left_obsB_blk = MT_good( [idx_all.st_left_obsB_bl] );
    
    %separate the above based on optimal vs average now
    MT_opt_right_obsA_tmp = MT_st_right_obsA_blk( [(opt_indc_st_right_obsA==1)] );
    MT_opt_left_obsA_tmp = MT_st_left_obsA_blk( [(opt_indc_st_left_obsA==1)] );
    
    MT_opt_right_obsB_tmp = MT_st_right_obsB_blk( [(opt_indc_st_right_obsB==1)] );
    MT_opt_left_obsB_tmp = MT_st_left_obsB_blk( [(opt_indc_st_left_obsB==1)] );
   
    %%%%MA
    MT_ma_right_obsA_tmp = MT_st_right_obsA_blk( [(opt_indc_st_right_obsA==0)] );
    MT_ma_left_obsA_tmp = MT_st_left_obsA_blk( [(opt_indc_st_left_obsA==0)] );
    
    MT_ma_right_obsB_tmp = MT_st_right_obsB_blk( [(opt_indc_st_right_obsB==0)] );
    MT_ma_left_obsB_tmp = MT_st_left_obsB_blk( [(opt_indc_st_left_obsB==0)] );
    
    figure(fn + 4); 
    
    %plot distribution of MTs
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_test, 'b', []);
%     if c3==1, title('All'); end
%     c3 = c3+1;
% %     xlim([0.4,2]);
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_st_obsA_blk, 'b', []);
%     if c3==2, title('STT obs A (blocked)'); end
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_st_obsB_blk, 'b', []);
%     if c3==3, title('STT obs B (blocked)'); end
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_dt_null_test, 'b', []);
%     if c3==4, title('DTT (null)'); end
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_dt_obsA, 'b', []);
%     if c3==5, title('DTT obs A (all)'); end
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_dt_obsB, 'b', []);
%     if c3==6, title('DTT obs B (all)'); end
%     c3 = c3+1;


    subplot(num_subjects, 6, c3); hold on;
    create_hist_init(MT_test, 'b', []);
    if c3==1, title('All'); end
    c3 = c3+1;
%     xlim([0.4,2]);
    
    subplot(num_subjects, 6, c3); hold on;
    create_hist_init_V2(MT_good(idx_all.st_left_obsA_bl), purple, []); hold on;
    create_hist_init_V2(MT_good(idx_all.st_right_obsA_bl), red_, []);
    if c3==2, title('STT obs A (blocked)'); end
    YL = ylim;
    mu = nanmean(MT_st_obsA_blk);
    h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
    %leg1=legend(h1);
    %set(leg1,'Location','Best');
    c3 = c3+1;
    
    subplot(num_subjects, 6, c3); hold on;
    create_hist_init_V2(MT_good(idx_all.st_left_obsB_bl), purple, []); hold on;
    create_hist_init_V2(MT_good(idx_all.st_right_obsB_bl), red_, []);
    if c3==3, title('STT obs B (blocked)'); end
    YL = ylim;
    mu = nanmean(MT_st_obsB_blk);
    h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
    %leg1=legend(h1);
    %set(leg1,'Location','Best');
    c3 = c3+1;
    
    subplot(num_subjects, 6, c3); hold on;
    create_hist_init(MT_dt_null_test, 'b', []);
    p_dtt_null_mt_all(kk) = prctile(MT_dt_null_test,95);
    if c3==4, title('DTT (null)'); end
    c3 = c3+1;
    
    subplot(num_subjects, 6, c3); hold on;
    create_hist_init_V2(MT_good([idx_all.dt_left_obsA_bl, idx_all.dt_left_obsA_nb,]), purple, []); hold on;
    create_hist_init_V2(MT_good([idx_all.dt_right_obsA_bl, idx_all.dt_right_obsA_nb]), red_, []);
    if c3==5, title('DTT obs A (all)'); end
    YL = ylim;
    mu = nanmean(MT_dt_obsA);
    h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
    %leg1=legend(h1);
    %set(leg1,'Location','Best');
    c3 = c3+1;
    
    subplot(num_subjects, 6, c3); hold on;
    create_hist_init_V2(MT_good([idx_all.dt_left_obsB_bl, idx_all.dt_left_obsB_nb,]), purple, []); hold on;
    create_hist_init_V2(MT_good([idx_all.dt_right_obsB_bl, idx_all.dt_right_obsB_nb]), red_, []);
     p_dtt_B_mt_all(kk) = prctile([MT_good([idx_all.dt_left_obsB_bl, idx_all.dt_left_obsB_nb,]); MT_good([idx_all.dt_right_obsB_bl, idx_all.dt_right_obsB_nb])],95);
    if c3==6, title('DTT obs B (all)'); end
    YL = ylim;
    mu = nanmean(MT_dt_obsB);
    h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
   % leg1=legend(h1);
    %set(leg1,'Location','Best');
    c3 = c3+1;
    
    
    
    figure(333); hold on;
    
%     MT_opt_right_obsA_tmp = MT_st_right_obsA_blk( [(opt_indc_st_right_obsA==1)] );
%     MT_opt_left_obsA_tmp = MT_st_left_obsA_blk( [(opt_indc_st_left_obsA==1)] );
%     
%     MT_opt_right_obsB_tmp = MT_st_right_obsB_blk( [(opt_indc_st_right_obsB==1)] );
%     MT_opt_left_obsB_tmp = MT_st_left_obsB_blk( [(opt_indc_st_left_obsB==1)] );

    
    
    %plot MT wrt trial number for STT blocked cases
    %label the tgt direction with color and label outward vs inward movements
    %with open vs closed circles
    
    %first plot the left target (purple) vs right (red) (make it 2 by 7) 
%     subplot(num_subjects,2,ss); hold on; grid minor;
%     if ~(isempty(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==1)) ))
%         plot(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==1)), MT_opt_right_obsA_tmp,'ro');
%     end
%     if ~(isempty(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==1)) ))
%         plot(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==1)), MT_opt_left_obsA_tmp,'o','color',purple);
%     end
%     
%    
%     if ~(isempty(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==0)) ))
%         plot(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==0)), MT_ma_right_obsA_tmp,'r.');
%     end
%     if ~(isempty(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==0)) ))
%         plot(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==0)), MT_ma_left_obsA_tmp,'.','color',purple);
%     end
%     
%     if ss==1, title('Obs A'); end
%     
%     ss = ss+1;
%     subplot(num_subjects,2,ss); hold on; grid minor;
%     if ~(isempty(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==1)) ))
%         plot(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==1)), MT_opt_right_obsB_tmp,'ro');
%     end
%     if ~(isempty(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==1)) ))
%         plot(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==1)), MT_opt_left_obsB_tmp,'o','color',purple);
%     end
%     
%     
%     if ~(isempty(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==0)) ))
%         plot(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==0)), MT_ma_right_obsB_tmp,'r.');
%     end
%     if ~(isempty(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==0)) ))
%         plot(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==0)), MT_ma_left_obsB_tmp,'.','color',purple);
%     end
% 
%     if ss==2, title('Obs B'); end
%     
%     ss = ss+1;


    subplot(2,num_subjects,ss); hold on; grid minor;
    if ~(isempty(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==1)) ))
        plot(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==1)), MT_opt_right_obsA_tmp,'ro');
    end
    if ~(isempty(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==1)) ))
        plot(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==1)), MT_opt_left_obsA_tmp,'o','color',purple);
    end
    
   
    if ~(isempty(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==0)) ))
        plot(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==0)), MT_ma_right_obsA_tmp,'r.');
    end
    if ~(isempty(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==0)) ))
        plot(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==0)), MT_ma_left_obsA_tmp,'.','color',purple);
    end
    
    %if ss==1, title('Obs A'); end
    
    
    subplot(2,num_subjects,ss2); hold on; grid minor;
    if ~(isempty(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==1)) ))
        plot(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==1)), MT_opt_right_obsB_tmp,'ro');
    end
    if ~(isempty(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==1)) ))
        plot(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==1)), MT_opt_left_obsB_tmp,'o','color',purple);
    end
    
    
    if ~(isempty(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==0)) ))
        plot(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==0)), MT_ma_right_obsB_tmp,'r.');
    end
    if ~(isempty(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==0)) ))
        plot(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==0)), MT_ma_left_obsB_tmp,'.','color',purple);
    end

   % if ss==2, title('Obs B'); end
    
    ss = ss+1;
    ss2 = ss2+1;
    
end


figure(fn+ 5);
subplot(2,2,1); hold on;
create_hist_init([MT_opt_obsA], grey, []);
title('obs A');
ylabel('optimal');
subplot(2,2,2); hold on;
create_hist_init([MT_ma_obsA], grey, []);
title('obs B');


subplot(2,2,3); hold on;
ylabel('averaging');
create_hist_init([MT_opt_obsB], grey, []);
subplot(2,2,4); hold on;
create_hist_init([MT_ma_obsB], grey, []);

tt1 = {['STTs for obs A'], ['STT for obs B'], ['DTT for obs A'], ['DTT for obs B'], ['Movement Times'], ['optimal vs averaging MT']};
for qz = 1:6
   figure(fn + (qz - 1) );
   suptitle(tt1{qz});
    %set(gca,'Fontsize',14);
end



%% Next we need to organize the data according to block
%Eventually, we will sample blocks for each subject in which optimal motions on STTs are made

tgt = info_all.tgt;
tpb = info_all.trials_per_block;

%fam (60 trials) -> fam (60 trials) -> ST obs (60 trials) -> ST obs (60 trials) -> DTT + STT baseline (60 trials) x 2
%-> obs 2 cm test (100 trials) * 4 -> obs 1 cm test (100 trials) * 4
block_ind = cell(1,length(tpb));
block_ind{1} = [1:tpb(1)];

%get indices of forward trials per block
for kk=2:length(tpb)
   block_ind{kk} =  [block_ind{kk-1}(end)+1:block_ind{kk-1}(end) + tpb(kk)]; 
end


%now we will pool the trial types indices in the sequence of the blocks
testb = 7; %this is the start of our test block
num_test = length(tpb) - testb +1;

% stt_obsA_bl.R = cell(num_subjects, num_test);
% stt_obsA_bl.L = stt_obsA_bl.R;
% 
% stt_obsB_bl.R = cell(num_subjects, num_test);
% stt_obsB_bl.L = stt_obsB_bl.R;
% 
% dtt_obsA = cell(num_subjects, num_test);
% %dtt_obsA.L = dtt_obsA_bl.R;
% 
% dtt_obsB = cell(num_subjects, num_test);
% %dtt_obsB.L = dtt_obsB_bl.R;

% dtt_obsA_right = nan(num_subjects, num_test);
% dtt_obsA_left = nan(num_subjects, num_test);
% dtt_obsB_right = nan(num_subjects, num_test);
% dtt_obsB_left = nan(num_subjects, num_test);

dtt_obsA_right = [];
dtt_obsA_left = [];
dtt_obsB_right = [];
dtt_obsB_left = [];

%%%%%
dt_obsA_right_idx = [idx_all.dt_right_obsA_bl, idx_all.dt_left_obsA_nb];
dt_obsA_left_idx = [idx_all.dt_left_obsA_bl, idx_all.dt_right_obsA_nb];

dt_obsB_right_idx = [idx_all.dt_right_obsB_bl, idx_all.dt_left_obsB_nb];
dt_obsB_left_idx = [idx_all.dt_left_obsB_bl, idx_all.dt_right_obsB_nb];

%indicator for which blocks to sample from which subjects
opt_block_sample = nan(num_subjects,num_test);


block_st_right_obsA_bl_idx = cell(1,num_test);
block_st_left_obsA_bl_idx = cell(1,num_test);

block_st_right_obsB_bl_idx = cell(1,num_test);
block_st_left_obsB_bl_idx = cell(1,num_test);

block_dt_obsA_right_idx = cell(1,num_test);
block_dt_obsA_left_idx = cell(1,num_test);
block_dt_obsB_right_idx = cell(1,num_test);
block_dt_obsB_left_idx = cell(1,num_test);


for jj=1:num_test
    
   current_block = jj + (testb) - 1;
   current_idx = block_ind{current_block};
   
   %get the indices we want
   %B needs to be the current_idx
    
   block_st_right_obsA_bl_idx{jj} = idx_all.st_right_obsA_bl(ismember(idx_all.st_right_obsA_bl, current_idx));
   block_st_left_obsA_bl_idx{jj} = idx_all.st_left_obsA_bl(ismember(idx_all.st_left_obsA_bl, current_idx));
   
   block_st_right_obsB_bl_idx{jj} = idx_all.st_right_obsB_bl(ismember(idx_all.st_right_obsB_bl, current_idx));
   block_st_left_obsB_bl_idx{jj} = idx_all.st_left_obsB_bl(ismember(idx_all.st_left_obsB_bl, current_idx));
   
   block_dt_obsA_right_idx{jj} = dt_obsA_right_idx(ismember( dt_obsA_right_idx, current_idx)); 
   block_dt_obsA_left_idx{jj} = dt_obsA_left_idx(ismember( dt_obsA_left_idx, current_idx)); 
   
   block_dt_obsB_right_idx{jj} = dt_obsB_right_idx(ismember( dt_obsB_right_idx, current_idx));
   block_dt_obsB_left_idx{jj} = dt_obsB_left_idx(ismember( dt_obsB_left_idx, current_idx)); 
   
   
   for kk=1:num_subjects
       %figure out if during certain blocks, a subject produced at least 70% optimal movements on STTs
       opt_idx = optimal_indicator_all(kk,:);
       
       block_opt_idx_obsA_right = opt_idx(block_st_right_obsA_bl_idx{jj});
       block_opt_idx_obsA_left = opt_idx(block_st_left_obsA_bl_idx{jj});
       
       block_opt_idx_obsB_right = opt_idx(block_st_right_obsB_bl_idx{jj});
       block_opt_idx_obsB_left = opt_idx(block_st_left_obsB_bl_idx{jj});
       
       
       %now determine if its a good block to sample
       tmp1 = sum(block_opt_idx_obsA_right) / length(block_opt_idx_obsA_right);
       tmp2 = sum(block_opt_idx_obsA_left) / length(block_opt_idx_obsA_left);
       tmp3 = sum(block_opt_idx_obsB_right) / length(block_opt_idx_obsB_right);
       tmp4 = sum(block_opt_idx_obsB_left) / length(block_opt_idx_obsB_left);
       
       tmp5 = [tmp1,tmp2,tmp3, tmp4];
       
       if any(tmp5> 0.7)
           opt_block_sample(kk,jj) = 1;
       else
           opt_block_sample(kk,jj) = 0;
       end         
   end
end

%figure(555); hold on;
%now lets get the data one want from the good blocks
for p = 1:num_subjects%6
    for q=1:num_test
        
        %sample the blocks we want
        if opt_block_sample(p,q)
            %get the angle!!! (and also trajectories???)
            
            for i=1:length(block_dt_obsA_right_idx{q})
                dtt_obsA_right = [dtt_obsA_right, initial_angle(p,block_dt_obsA_right_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
                
                 %%%check some random trials
%                  if initial_angle(p,block_dt_obsA_right_idx{q}(i)) <= -60
%                      keyboard;
%                  end
                
            end
            
            for i=1:length(block_dt_obsA_left_idx{q})
                dtt_obsA_left = [dtt_obsA_left, initial_angle(p,block_dt_obsA_left_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
            end
            
            for i=1:length(block_dt_obsB_right_idx{q})
                dtt_obsB_right = [dtt_obsB_right, initial_angle(p,block_dt_obsB_right_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
            end
            
            for i=1:length(block_dt_obsB_left_idx{q})
                dtt_obsB_left = [dtt_obsB_left, initial_angle(p,block_dt_obsB_left_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
            end
            
%             if p==6,
%                 figure(555); hold on;
%                 for zz=1:length(block_dt_obsA_left_idx{q})
%                     plot(initial_mov{p,block_dt_obsA_left_idx{q}(zz)}(:,1), initial_mov{p,block_dt_obsA_left_idx{q}(zz)}(:,2) );
%                 end
%             end
            
            %if p==6 & q==3, keyboard; end
            
            
        end
    end
end



figure(fn+ 6);
subplot(2,2,1); hold on;
create_hist_init(dtt_obsA_right, grey, []);
YL = ylim;
mu = nanmean(dtt_obsA_right);
h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
%leg1=legend(h1);
%set(leg1,'Location','Best');
title('Right obstacle');
ylabel('obs A');
subplot(2,2,2); hold on;
create_hist_init(dtt_obsA_left, grey, []);
title('Left Obstacle');
YL = ylim;
mu = nanmean(dtt_obsA_left);
h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% leg1=legend(h1);
% set(leg1,'Location','Best');


subplot(2,2,3); hold on;
ylabel('obs B');
create_hist_init(dtt_obsB_right, grey, []);
YL = ylim;
mu = nanmean(dtt_obsB_right);
h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% leg1=legend(h1);
% set(leg1,'Location','Best');

subplot(2,2,4); hold on;
create_hist_init(dtt_obsB_left, grey, []);
YL = ylim;
mu = nanmean(dtt_obsB_left);
h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% leg1=legend(h1);
% set(leg1,'Location','Best');

%%

%pool angles across all participants from DTT trials and make a histogram
dt_iang_pooled = [];
for i=1:num_subjects
%dt_iang_pooled = [right_obsB_dtt_iang(i,:), left_obsB_dtt_iang(i,:) * -1, dt_iang_pooled];
dt_iang_pooled = [right_obsB_dtt_iang(i,:), dt_iang_pooled];
end

figure; hold on; title('distribution of angles across all subjects');
create_hist_init(dt_iang_pooled, grey, 40);
YL = ylim;
mu = nanmedian(dt_iang_pooled);
h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% leg1=legend(h1);
% set(leg1,'Location','Best');

%% Organize the ST data
obs_fields = fieldnames(inside_indc);
%fields are (in order): obsA_L, obsA_R, obsB_L, obsB_R

%will also optionally account for baseline subtraction
st_block_dat = struct( obs_fields{1}, iang_st_left_obsA_bl-30 - (ones(size(iang_st_left_obsA_bl))' * diag(bias_st_left) * st_bias_flag)', ...
    obs_fields{2}, iang_st_right_obsA_bl+30 - (ones(size(iang_st_right_obsA_bl))' * diag(bias_st_right) * st_bias_flag)',...
    obs_fields{3}, iang_st_left_obsB_bl-30 - (ones(size(iang_st_left_obsB_bl))' * diag(bias_st_left) * st_bias_flag)',...
    obs_fields{4}, iang_st_right_obsB_bl+30 - (ones(size(iang_st_right_obsB_bl))' * diag(bias_st_right) * st_bias_flag)' );


%% Optionally remove inside ST trials before the prediction analysis is done 
%only remove trials for subjects who didnt exhibit inside behavior on the majority of trials

inside_filter = 1;

if inside_filter
    fn = fn+50;
    
    for qo=1:length(obs_fields)
        
        cfld = obs_fields{qo};
        
        bad_mat = double(~(inside_indc.(cfld)'));
        bad_mat(bad_mat==0) = nan;
        
        st_block_dat.(cfld) = st_block_dat.(cfld) .* bad_mat;
    end
end

%%
%close all;

%make a 2 x2 scatter plot of the average of the median angles for STT
%left and right targets for each subject (motor averaging prediction)
%vs the observed angle on DTTs
Lw1 = 1.5;
red_=[1 0 0];
purple=[148,0,211]/256;
col1 = red_;
col2 = purple;
xq = [-50:1:50];
yq =xq/2;
sf = 1;

figure(fn + 10); hold on;
%first is obs A Right target

avg_pred_obsA_right = ( (nanmedian(iang_st_left_obsA_nb,2) - bias_st_left *st_bias_flag) + ...
    (nanmedian(iang_st_right_obsA_bl,2) -bias_st_right *st_bias_flag) ) /2;
avg_obsA_right =  (nanmedian(iang_st_right_obsA_bl,2) -bias_st_right *st_bias_flag)+30;
dtt_obsA_right_all = nanmedian(right_obsA_dtt_iang,2);

avg_pred_obsA_left = ( (nanmedian(iang_st_right_obsA_nb,2) - bias_st_right *st_bias_flag) + ...
    (nanmedian(iang_st_left_obsA_bl,2) - bias_st_left *st_bias_flag) ) /2;
avg_obsA_left = nanmedian(iang_st_left_obsA_bl,2) - bias_st_left *st_bias_flag - 30;
dtt_obsA_left_all = nanmedian(left_obsA_dtt_iang,2);

avg_pred_obsB_right = ( (nanmedian(iang_st_left_obsB_nb,2) -bias_st_left *st_bias_flag) + ...
    (nanmedian(iang_st_right_obsB_bl,2)-bias_st_right *st_bias_flag) ) /2;
avg_obsB_right =  (nanmedian(iang_st_right_obsB_bl,2)-bias_st_right *st_bias_flag) +30;
dtt_obsB_right_all = nanmedian(right_obsB_dtt_iang,2);

avg_pred_obsB_left = ( (nanmedian(iang_st_right_obsB_nb,2) -bias_st_right *st_bias_flag) + ...
    (nanmedian(iang_st_left_obsB_bl,2)-bias_st_left *st_bias_flag) ) /2;
avg_obsB_left = (nanmedian(iang_st_left_obsB_bl,2)-bias_st_left *st_bias_flag)-30 ;
dtt_obsB_left_all = nanmedian(left_obsB_dtt_iang,2);

% phi1 = 12.5288;
% phi2 = 6.3402;
phiA = atand(2/9); %how much the obstacle sticks out, and where it is frist seen along its axis
phiB = atand(1/9); %for obstacle B

if sf==0
    subplot(2,2,1); hold on;
    
    plot(nanmean(avg_pred_obsA_right), nanmean(dtt_obsA_right_all),'p','markersize',12,'linewidth',Lw1,'color',col1);
    
    plot(xq,yq,'--','color','k');
    %plot(-xq,yq,'--','color',grey);
    
    plot(avg_pred_obsA_right, dtt_obsA_right_all, 'o', 'linewidth' ,Lw1,'color',col1);
    xlabel('observed angle from STT');
    ylabel('observed angle on DTT');
    title('obs A, right side');
    plot_ellipse_10_30(avg_pred_obsA_right,dtt_obsA_right_all,col1,Lw1);
    
    %now do obs A left
    subplot(2,2,2); hold on;
    plot(xq,yq,'--','color','k');
    %plot(-xq,yq,'--','color',grey);
    plot(avg_pred_obsA_left, dtt_obsA_left_all, 'o', 'linewidth' ,Lw1,'color',col1);
    plot(nanmean(avg_pred_obsA_left), nanmean(dtt_obsA_left_all),'p','markersize',12,'linewidth',Lw1,'color',col1);
    xlabel('observed angle from STT');
    ylabel('observed angle on DTT');
    title('obs A, left side');
    plot_ellipse_10_30(avg_pred_obsA_left,dtt_obsA_left_all,col1,Lw1);
    
    
    %next is obs B right
    subplot(2,2,3); hold on;
    
    plot(nanmean(avg_pred_obsB_right), nanmean(dtt_obsB_right_all),'p','markersize',12,'linewidth',Lw1,'color',col2);
    plot(xq,yq,'--','color','k');
    %plot(-xq,yq,'--','color',grey);
    plot(avg_pred_obsB_right, dtt_obsB_right_all, 'o', 'linewidth' ,Lw1,'color',col2);
    plot_ellipse_10_30(avg_pred_obsB_right,dtt_obsB_right_all,col2,Lw1);
    xlabel('observed angle from STT');
    ylabel('observed angle on DTT');
    title('obs B, right side');
    
    
    %next is obs B left
    subplot(2,2,4); hold on;
    
    plot(nanmean(avg_pred_obsB_left), nanmean(dtt_obsB_left_all),'p','markersize',12,'linewidth',Lw1,'color',col2);
    plot(xq,yq,'--','color','k');
    %plot(-xq,yq,'--','color',grey);
    plot(avg_pred_obsB_left, dtt_obsB_left_all, 'o', 'linewidth' ,Lw1,'color',col2);
    plot_ellipse_10_30(avg_pred_obsB_left,dtt_obsB_left_all,col2,Lw1);
    xlabel('observed angle from STT');
    ylabel('observed angle on DTT');
    title('obs B, left side');
    
    for k=1:4, subplot(2,2,k); qk=30; xlim(qk*[-1, 1]); ylim(qk*[-1,1]); plot([0 0],qk*[-1,1],'k');plot(qk*[-1,1],[0 0],'k'); end; shg
    
else
    
    test_shade=[230,230,230]/256;
    
    %do the predictions
    title('Obs A Predictions');
    %start with obs A, left
    display_optimal_prediction_V1( st_block_dat.obsA_L, left_obsA_dtt_iang, -1, phiA, []);
    
    %obs A, right
    display_optimal_prediction_V1( st_block_dat.obsA_R, right_obsA_dtt_iang, 1, phiA, []);
    
   %now do the same for obstacle B on a new figure
    figure(fn + 11); hold on; 
    title('Obs B Predictions');
    
    %obs B, left
    display_optimal_prediction_V1( st_block_dat.obsB_L, left_obsB_dtt_iang, -1, phiB, []);
    
    %obs B, right
    display_optimal_prediction_V1( st_block_dat.obsB_R, right_obsB_dtt_iang, 1, phiB, []);
    
    
end

%% make a time series plot of movement times, mark rewarded vs unrewarded
%make one for test and 2 others for baseline STT w/ obstacles (mark blocked
%obstacles) and DTT blocks as well
% qr = 1;
% 
% for jk=1:num_subjects
%     
%     MT_label_sub = MT_labels_sub_all(jk,:)'; % 1 means it was too slow
%     MT = dat_all.MT_sub_all(jk,:);
%     
%     TF = find(MT <= .225); %too fast
%     MT_label_sub(TF) = -1;
%     
%     %ignore the glass breaking stuff since it might have been wrong
%     sub_bad_trial = squeeze(info_all.Bad_labels_sub_all(kk,:,:));
%     sub_good_trial = ~sub_bad_trial;
%     
%     good_label= prod(sub_good_trial,2);
%     MT_good = MT'.*good_label;
%     
%     %first the test periods
%     figure(111);
%     %2 columns for each obstacle type, for now look at this during test period
%     subplot(num_subjects, 2,qr+(jk-1)); hold on;
%     test1_idx = exp_seq{7};
%     MT_test = MT_good(test1_idx);
%     
%     display_MT_rwd(MT_test,test1_idx,MT_label_sub(test1_idx));
%     qr = qr+1;
%     
%     
%     subplot(num_subjects, 2,qr+(jk-1)); hold on;
%     display_MT_rwd(MT_good(exp_seq{8}),exp_seq{8},MT_label_sub(exp_seq{8}));
%     %qr = qr+1;
% end

%%

%close all;

%for instances where a subject did an inside movement (during test), lets look at the
%previous N movements of that type to see if they're rewarded (STT blocked trials)

num_prev_rwd = 3; %look at the three previous rewarded movements

stt_block_idx.obsA_L = idx_all.st_left_obsA_bl;
stt_block_idx.obsA_R = idx_all.st_right_obsA_bl;
stt_block_idx.obsB_L = idx_all.st_left_obsB_bl;
stt_block_idx.obsB_R = idx_all.st_right_obsB_bl;

%combine and save all the blocked indices together (used for later)
st_block_idx_comb = [idx_all.st_left_obsA_bl, idx_all.st_right_obsA_bl, idx_all.st_left_obsB_bl, idx_all.st_right_obsB_bl];

prev_rwd = nan(num_prev_rwd+1, num_subjects, length(obs_fields));

for k=1:num_subjects
   
  for j=1:length(obs_fields)
      
      cf = obs_fields{j};
      
      cidx = inside_indc.(cf)(:,k);
      
      if any(cidx) %there was an outside trial
          
          inside_idx_tmp = find(cidx==1, 1, 'last'); %for now, just focus on finding one bad trial per subject
          
          %find the previous three trials
          try
              rwd_idx_tmp = stt_block_idx.(cf)(max( 1, inside_idx_tmp - num_prev_rwd): inside_idx_tmp);
          catch
              keyboard; end
          
          prev_rwd([1:length(rwd_idx_tmp)],k, j) = sub_rwd_all(rwd_idx_tmp, k);
          
          
      end
      
      
  end
    
    
end

% figure; hold on;
% title('obs A rwd for all subjects, left');
% imagesc(prev_rwd(:,:,1));
% myColorMap = jet;
% myColorMap(1,:) = 1;
% myColorMap(end,:) = 0;
% colormap(myColorMap);

figure; title('Obs A rwd for all subjects'); 
imshow(prev_rwd(:,:,1), 'InitialMagnification', 'fit');

% figure; hold on;
% title('obs A rwd for all subjects, right');
% imagesc(prev_rwd(:,:,2));
% myColorMap = jet;
% myColorMap(1,:) = 1;
% myColorMap(end,:) = 0;
% colormap(myColorMap);


%% also look at how these movements are rewarded in aggregate, and then look at rwd rates for stt blocked trials vs all the other ones (combined)!!

st_block_test_rwd_all = struct(obs_fields{1}, [], obs_fields{2}, [], obs_fields{3}, [], obs_fields{4}, [] );
st_block_test_rwd_all.comb = [];
%sub_rwd_all_good = double(~sub_rwd_all); %now a 1 is good

for qz =1:num_subjects
    for w1=1:length(obs_fields)
       rwd_trials_tmp = rwd_label(qz, stt_block_idx.(obs_fields{w1}) );
       st_block_test_rwd_all.(obs_fields{w1})(qz) = sum(rwd_trials_tmp) / length(rwd_trials_tmp) * 100;
    end
end

%look at blocked trials in aggregate
for qz2=1:num_subjects
    st_block_test_rwd_all.comb(qz2) = sum(rwd_label(qz2, st_block_idx_comb)) / length(st_block_idx_comb) *100;
end


%look at all stt trials during the test phase combined
st_all_score = nan(num_subjects, 1);
dt_all_score = nan(num_subjects, 1);
for qz3 = 1:num_subjects
    st_all_score(qz3) = sum(rwd_label(qz3, idx_all.st_test_all)) / length(idx_all.st_test_all) * 100;
    dt_all_score(qz3) = sum(rwd_label(qz3, idx_all.dt_test_all)) / length(idx_all.dt_test_all) * 100;
end

%now see how the scores fare when we remove the blocked obs trials
st_idx_all_tmp = find(ismember( idx_all.st_test_all, st_block_idx_comb));
st_idx_all_no_obs = idx_all.st_test_all;
st_idx_all_no_obs(st_idx_all_tmp) = [];

st_no_obs_score = nan(num_subjects, 1);

for qz4 = 1:num_subjects
    st_no_obs_score(qz4) = sum(rwd_label(qz4, st_idx_all_no_obs)) / length(st_idx_all_no_obs) * 100;
end

%plot the scores per block
figure; hold on;
for qz5 = 1:num_subjects
   plot( info_all.Score_sub_all(qz5,:));
end

xlabel('exp block');
ylabel('Score');

%could probably also make a time series of the feedback given (check it
%early on especially...)

%% early vs late analysis

% do the predictions we have above, but this time only use early and late samples

early_frac = 0.25; %fraction of trials to use from the beginning
late_frac = 0.25; %same for the late period

%also divide the data into quarters to see if that has an effect
idx_st_all = reshape([1:length(iang_st_left_obsA_bl)], length(iang_st_left_obsA_bl)/4, 4);

early_period_obsA_idx = idx_st_all(:,3);
%early_period_obsB_idx = early_idx_st_all(:,1);
late_period_obsA_idx = idx_st_all(:,4);

% early_period_obsA_idx = [1: floor(length(iang_st_left_obsA_bl) * early_frac)];
 early_period_obsB_idx = [1: floor(length(iang_st_left_obsB_bl) * early_frac)]; %last quarter for these are only 2 trials!

%late_period_obsA_idx = [length(iang_st_left_obsA_bl) - floor(length(iang_st_left_obsA_bl) * late_frac) : length(iang_st_left_obsA_bl)];
late_period_obsB_idx = [length(iang_st_left_obsB_bl) - floor(length(iang_st_left_obsB_bl) * late_frac) : length(iang_st_left_obsB_bl)];

early_dat.st_left_A = st_block_dat.obsA_L(:,early_period_obsA_idx);
early_dat.st_right_A = st_block_dat.obsA_R(:,early_period_obsA_idx);
early_dat.st_left_B = st_block_dat.obsB_L(:,early_period_obsB_idx);
early_dat.st_right_B = st_block_dat.obsB_R(:,early_period_obsB_idx);

early_dat.dt_left_A = left_obsA_dtt_iang(:, early_period_obsA_idx); %since there are more dt trials, this is not actually using the last quarter!
early_dat.dt_right_A = right_obsA_dtt_iang(:, early_period_obsA_idx);
early_dat.dt_left_B = left_obsB_dtt_iang(:, early_period_obsB_idx);
early_dat.dt_right_B = right_obsB_dtt_iang(:, early_period_obsB_idx);

late_dat.st_left_A = st_block_dat.obsA_L(:,late_period_obsA_idx);
late_dat.st_right_A = st_block_dat.obsA_R(:,late_period_obsA_idx);
late_dat.st_left_B = st_block_dat.obsB_L(:,late_period_obsB_idx);
late_dat.st_right_B = st_block_dat.obsB_R(:,late_period_obsB_idx);

late_dat.dt_left_A = left_obsA_dtt_iang(:, late_period_obsA_idx);
late_dat.dt_right_A = right_obsA_dtt_iang(:, late_period_obsA_idx);
late_dat.dt_left_B = left_obsB_dtt_iang(:, late_period_obsB_idx);
late_dat.dt_right_B = right_obsB_dtt_iang(:, late_period_obsB_idx);


%%%%%%%%%%% plot them
figure; hold on; title('Early vs Late: Obs A');

%do it for obs A first
display_optimal_prediction_EL( early_dat.st_left_A, early_dat.dt_left_A, late_dat.st_left_A, late_dat.dt_left_A, -1, phiA, [] );
display_optimal_prediction_EL( early_dat.st_right_A, early_dat.dt_right_A, late_dat.st_right_A, late_dat.dt_right_A, 1, phiA, [] );

%now for obstacle B
figure; hold on; title('Early vs Late: Obs B');

display_optimal_prediction_EL( early_dat.st_left_B, early_dat.dt_left_B, late_dat.st_left_B, late_dat.dt_left_B, -1, phiB, [] );
display_optimal_prediction_EL( early_dat.st_right_B, early_dat.dt_right_B, late_dat.st_right_B, late_dat.dt_right_B, 1, phiB, [] );


%% 













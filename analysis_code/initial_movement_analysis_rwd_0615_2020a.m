%additional analysis for the obstacle stuff
%run the other analysis before this one!!
close all;
load('group_analysis.mat');
clearvars -except Cdr dat_all exp_seq Ideal_dir_all idx_all info_all MT_labels_sub_all...
    num_subjects optimal_indicator_all Pos_sub_all S Sub_id_all xTraj_sub_all yTraj_sub_all...
    rwd_label

tgt = info_all.tgt;

red_=[1 0 0];
purple=[128,0,128]/256;
violet = [128, 130, 238] / 256;
grey=[0.5 0.5 0.5];
orange=[255,140,0]/256;
black_=[0,0,0];
dblue = [0,0,255] / 256;
cyan = [0,255,255] / 256;
light_blue = [0,0,255]/ 255;

%define the obstacle location (in normalized coorindates)
obs_locx = 0.5 * 20 * cosd(60);
obs_locy = 0.5 * 20 * sind(60);

%calculate the angle for all the cases at a certain point into the movement

initial_marker = 10; %cm into the movement

%do we want to subtract by individual subject biases?
st_bias_flag = 0;
dt_bias_flag = 0; 

inside_filter = 1; %option to remove trials where subjects move inside relative to the obstacle on STT blocked trials

iqr_filter_flag = 0; %flag to turn on filtering of data based on iqr per subject, turn off if already ran in group analysis script

%assuming that all subjects completed the same number of trials
num_trials = length(Pos_sub_all(1,:));

%calculate the angle and save it
initial_angle = nan(num_subjects, num_trials);
%v1 = [0,initial_marker]; %vector we will calculate the angle wrt

initial_mov = cell(num_subjects, num_trials);

for qq = 1:num_subjects
    for kn = 1:num_trials
        
        x_data = (xTraj_sub_all{qq,kn});
        y_data = yTraj_sub_all{qq,kn};
        
        current_pos = Pos_sub_all{qq,kn};
        %current_pos = y_data-y_data(1);
        
        %determine if its a 1-tgt or 2-tgt trial (1 or 3 on column 11)
        %if its a 1-target, if x position of target (column 3) is bigger,
        %less, or equal to 950, then its a right, left, and straight target
        %respectively
        if tgt(kn,11)==3
            %initial_marker = 10; %use 10 cm for 2-tgt trials
            initial_marker = obs_locy - 2*sind(60);
        else
            initial_marker = obs_locy - 2*sind(60); %get the y position of the end of the obstacle
%             if tgt(kn,3) > 950 %right target
%                 initial_marker = obs_locy;
%             elseif tgt(kn,3) < 950 %left
%                 isLeft(kk)=1;
%             elseif tgt(kn,3) == 950 %straight
%                 isStraight(kk)=1;
%             end    
        end
        
        %optionally, we will want to calculate the angle at axis of the obstacle
        %this will depend on the target location for the given trial type
        
        initial_idx = find(current_pos >= initial_marker,1,'first')-2; %take an extra sample back just to be safe
        %now use this to get the vector for angle calculations
        
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

% xtest_null_dt = xTraj_sub_all(:, idx_all.dt_null_test );
% ytest_null_dt = yTraj_sub_all(:, idx_all.dt_null_test );

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

% xtest_null_dt_tmp = cellfun(@(x) x(1:max_samples), xtest_null_dt,'UniformOutput',false);
% ytest_null_dt_tmp = cellfun(@(x) x(1:max_samples), ytest_null_dt,'UniformOutput',false);

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
%iang_dt_null_test = nan(num_subjects, 1);

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
    
    %xtest_dt_null_avg(jk,:) =  nanmedian( cell2mat(xtest_null_dt_tmp(jk,:)), 2);
    %ytest_dt_null_avg(jk,:) =  nanmedian( cell2mat(ytest_null_dt_tmp(jk,:)), 2);
   
   %calculate the angle based on the average trajectory
   P1 = cumsum(sqrt(diff(xbln_st_left_avg(jk,:)).^2+diff(ybln_st_left_avg(jk,:)).^2));
   P2 = cumsum(sqrt(diff(xbln_st_right_avg(jk,:)).^2+diff(ybln_st_right_avg(jk,:)).^2));
   P3 = cumsum(sqrt(diff(xbln_st_strt_avg(jk,:)).^2+diff(ybln_st_strt_avg(jk,:)).^2));
   P4 = cumsum(sqrt(diff(xbln_dt_avg(jk,:)).^2+diff(ybln_dt_avg(jk,:)).^2));
   %P5 = cumsum(sqrt(diff(xtest_dt_null_avg(jk,:)).^2+diff(ytest_dt_null_avg(jk,:)).^2));
   
   %NOTE: change above to ypos!!!
   
   idx1 = find(P1 >=initial_marker,1,'first');
   idx2 = find(P2 >=initial_marker,1,'first');
   idx3 = find(P3 >=initial_marker,1,'first');
   idx4 = find(P4 >=initial_marker,1,'first');
   %idx5 = find(P5 >=initial_marker,1,'first');
   
   v1_st_left = [xbln_st_left_avg(jk,1) - S(1), initial_marker];
   v1_st_right = [xbln_st_right_avg(jk,1) - S(1), initial_marker];
   v1_st_strt = [xbln_st_strt_avg(jk,1) - S(1),initial_marker];
   v1_dt = [xbln_dt_avg(jk,1) - S(1),initial_marker];
   
   v2_st_left = [ xbln_st_left_avg(jk,idx1) - S(1), ybln_st_left_avg(jk,idx1) - S(2) ];
   v2_st_right = [ xbln_st_right_avg(jk,idx2) - S(1), ybln_st_right_avg(jk,idx2) - S(2) ];
   v2_st_strt = [ xbln_st_strt_avg(jk,idx3) - S(1), ybln_st_strt_avg(jk,idx3) - S(2) ];
   v2_dt = [ xbln_dt_avg(jk,idx4) - S(1), ybln_dt_avg(jk,idx4) - S(2) ];
   %v2_dt_null_test = [ xtest_dt_null_avg(jk,idx5) - S(1), ytest_dt_null_avg(jk,idx5) - S(2) ];
   
   %get the angle
   iang_st_left_bln(jk) = atan2d(det([v1_st_left;v2_st_left]),dot(v1_st_left,v2_st_left));
   iang_st_right_bln(jk) = atan2d(det([v1_st_right;v2_st_right]),dot(v1_st_right,v2_st_right));
   iang_st_strt_bln(jk) = atan2d(det([v1_st_strt;v2_st_strt]),dot(v1_st_strt,v2_st_strt));
   iang_dt_bln(jk) = atan2d(det([v1_dt; v2_dt]),dot(v1_dt, v2_dt));
   %iang_dt_null_test(jk) = atan2d(det([v1;v2_dt_null_test]),dot(v1,v2_dt_null_test));
   
end
        
%now we can calculate the bias...
%using the average trajectory might not be as robust as I think (without interpolating!)
bias_st_left = (iang_st_left_bln - 30) * 1;
bias_st_right = (iang_st_right_bln + 30) * 1;
bias_st_strt = iang_st_strt_bln + 0;
bias_dt_null_bln = iang_dt_bln + 0;

bias_dt_null = bias_dt_null_bln;


%% now index these for all our cases

iang_st_left_obsA_bl = nan(num_subjects, size(idx_all.st_left_obsA_bl,2) );
iang_st_right_obsA_bl = nan(num_subjects, size(idx_all.st_right_obsA_bl,2) );
iang_st_left_obsB_bl = nan(num_subjects, size(idx_all.st_left_obsB_bl,2) );
iang_st_right_obsB_bl = nan(num_subjects, size(idx_all.st_right_obsB_bl,2) );

iang_st_left_obsA_nb = nan(num_subjects, size(idx_all.st_left_obsA_nb,2) );
iang_st_right_obsA_nb = nan(num_subjects, size(idx_all.st_right_obsA_nb,2) );
iang_st_left_obsB_nb = nan(num_subjects, size(idx_all.st_left_obsB_nb,2) );
iang_st_right_obsB_nb = nan(num_subjects, size(idx_all.st_right_obsB_nb,2) );

iang_dt_left_obsA_bl = nan(num_subjects, size(idx_all.dt_left_obsA_bl,2) );
iang_dt_right_obsA_bl = nan(num_subjects, size(idx_all.dt_right_obsA_bl,2) );
iang_dt_left_obsB_bl = nan(num_subjects, size(idx_all.dt_left_obsB_bl,2) );
iang_dt_right_obsB_bl = nan(num_subjects, size(idx_all.dt_right_obsB_bl,2) );

iang_dt_left_obsA_nb = nan(num_subjects, size(idx_all.dt_left_obsA_nb,2) );
iang_dt_right_obsA_nb = nan(num_subjects, size(idx_all.dt_right_obsA_nb,2) );
iang_dt_left_obsB_nb = nan(num_subjects, size(idx_all.dt_left_obsB_nb,2) );
iang_dt_right_obsB_nb = nan(num_subjects, size(idx_all.dt_right_obsB_nb,2) );

iang_strt_obsA_left = nan(num_subjects, size(idx_all.strt_obsA_left,2) );
iang_strt_obsA_right = nan(num_subjects, size(idx_all.strt_obsA_right,2) );
iang_strt_obsB_left = nan(num_subjects, size(idx_all.strt_obsB_left,2) );
iang_strt_obsB_right = nan(num_subjects, size(idx_all.strt_obsB_right,2) );

iang_dt_bln = nan(num_subjects, size(idx_all.dt_bln,2));
%idx_all.st_left_bln
iang_st_left_bln = nan(num_subjects, size(idx_all.st_left_bln,2));
iang_st_right_bln = nan(num_subjects, size(idx_all.st_right_bln,2));
iang_st_strt_bln = nan(num_subjects, size(idx_all.st_strt_bln,2));

% first for the angles
for qs = 1:num_subjects
    iang_st_left_obsA_bl(qs,:) = initial_angle(qs,idx_all.st_left_obsA_bl(qs,:));
    iang_st_right_obsA_bl(qs,:) = initial_angle(qs,idx_all.st_right_obsA_bl(qs,:));
    iang_st_left_obsB_bl(qs,:) = initial_angle(qs,idx_all.st_left_obsB_bl(qs,:));
    iang_st_right_obsB_bl(qs,:) = initial_angle(qs,idx_all.st_right_obsB_bl(qs,:));
    
    iang_st_left_obsA_nb(qs,:) = initial_angle(qs,idx_all.st_left_obsA_nb(qs,:));
    iang_st_right_obsA_nb(qs,:) = initial_angle(qs,idx_all.st_right_obsA_nb(qs,:));
    iang_st_left_obsB_nb(qs,:) = initial_angle(qs,idx_all.st_left_obsB_nb(qs,:));
    iang_st_right_obsB_nb(qs,:) = initial_angle(qs,idx_all.st_right_obsB_nb(qs,:));
    
    iang_dt_left_obsA_bl(qs,:) = initial_angle(qs,idx_all.dt_left_obsA_bl(qs,:));
    iang_dt_right_obsA_bl(qs,:) = initial_angle(qs,idx_all.dt_right_obsA_bl(qs,:));
    iang_dt_left_obsB_bl(qs,:) = initial_angle(qs,idx_all.dt_left_obsB_bl(qs,:));
    iang_dt_right_obsB_bl(qs,:) = initial_angle(qs,idx_all.dt_right_obsB_bl(qs,:));
    
    iang_dt_left_obsA_nb(qs,:) = initial_angle(qs,idx_all.dt_left_obsA_nb(qs,:));
    iang_dt_right_obsA_nb(qs,:) = initial_angle(qs,idx_all.dt_right_obsA_nb(qs,:));
    iang_dt_left_obsB_nb(qs,:)= initial_angle(qs,idx_all.dt_left_obsB_nb(qs,:));
    iang_dt_right_obsB_nb(qs,:) = initial_angle(qs,idx_all.dt_right_obsB_nb(qs,:));
    
    iang_strt_obsA_left(qs,:) = initial_angle(qs,idx_all.strt_obsA_left(qs,:));
    iang_strt_obsA_right(qs,:) = initial_angle(qs,idx_all.strt_obsA_right(qs,:));
    iang_strt_obsB_left(qs,:) = initial_angle(qs,idx_all.strt_obsB_left(qs,:));
    iang_strt_obsB_right(qs,:) = initial_angle(qs,idx_all.strt_obsB_right(qs,:));
    
    iang_dt_bln(qs,:) = initial_angle(qs, idx_all.dt_bln(qs,:));
    iang_st_left_bln(qs,:) = initial_angle(qs,idx_all.st_left_bln(qs,:));
    iang_st_right_bln(qs,:) = initial_angle(qs,idx_all.st_right_bln(qs,:));
    iang_st_strt_bln(qs,:) = initial_angle(qs,idx_all.st_strt_bln(qs,:));
end

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
    
    opt_indc_st_left_obsA = opt_idx(idx_all.st_left_obsA_bl(kk,:)) ;
    opt_indc_st_right_obsA = opt_idx(idx_all.st_right_obsA_bl(kk,:));
    opt_indc_st_left_obsB = opt_idx(idx_all.st_left_obsB_bl(kk,:));
    opt_indc_st_right_obsB = opt_idx(idx_all.st_right_obsB_bl(kk,:));
    
    opt_st_left_obsA = initial_angle(kk,idx_all.st_left_obsA_bl(kk,opt_indc_st_left_obsA==1) ) - (bias_st_left(kk) * st_bias_flag);
    opt_st_right_obsA = initial_angle(kk,idx_all.st_right_obsA_bl(kk,opt_indc_st_right_obsA==1) ) - (bias_st_right(kk) * st_bias_flag);
    opt_st_left_obsB = initial_angle(kk,idx_all.st_left_obsB_bl(kk,opt_indc_st_left_obsB==1) ) - (bias_st_left(kk) * st_bias_flag);
    opt_st_right_obsB = initial_angle(kk,idx_all.st_right_obsB_bl(kk,opt_indc_st_right_obsB==1) ) - (bias_st_right(kk) * st_bias_flag);
    
    ma_st_left_obsA = initial_angle(kk,idx_all.st_left_obsA_bl(kk,opt_indc_st_left_obsA==0) ) - (bias_st_left(kk) * st_bias_flag);
    ma_st_right_obsA = initial_angle(kk,idx_all.st_right_obsA_bl(kk,opt_indc_st_right_obsA==0) ) - (bias_st_right(kk) * st_bias_flag);
    ma_st_left_obsB = initial_angle(kk,idx_all.st_left_obsB_bl(kk,opt_indc_st_left_obsB==0) ) - (bias_st_left(kk) * st_bias_flag);
    ma_st_right_obsB = initial_angle(kk,idx_all.st_right_obsB_bl(kk,opt_indc_st_right_obsB==0) ) - (bias_st_right(kk) * st_bias_flag);
    
    %save the inside vs outside indices for STT trials, when people when
    %outside, we might not want to count this in our analysis
    inside_indc.obsA_L(:,kk) = ~opt_indc_st_left_obsA; % 1 means they didnt do the optimal (intended) movement
    inside_indc.obsA_R(:,kk) = ~opt_indc_st_right_obsA;
    inside_indc.obsB_L(:,kk) = ~opt_indc_st_left_obsB;
    inside_indc.obsB_R(:,kk) = ~opt_indc_st_right_obsB;
    
    
%     if kk==5, keyboard; end
    
%     figure(fn);
%     subplot(num_subjects,2,kk+c1); hold on;
%     create_hist_init(opt_st_left_obsA ,violet, 30);
%     create_hist_init(ma_st_left_obsA ,purple, 30);
%     if (kk+c1) == 1, title('Left target'); end
%     
%     subplot(num_subjects,2,kk+c2); hold on;
%     create_hist_init(opt_st_right_obsA ,orange, 30);
%     create_hist_init(ma_st_right_obsA ,red_, 30);
%     if (kk+c2) == 2, title('Right target'); end
%     
%     %repeat for obstacle B
%     figure(fn+1);
%     subplot(num_subjects,2,kk+c1); hold on;
%     create_hist_init(opt_st_left_obsB ,violet, 30);
%     create_hist_init(ma_st_left_obsB ,purple, 30);
%     if (kk+c1) == 1, title('Left target'); end
%     
%     
%     subplot(num_subjects,2,kk+c2); hold on;
%     create_hist_init(opt_st_right_obsB ,orange, 30);
%     create_hist_init(ma_st_right_obsB ,red_, 30);
%     if (kk+c2) == 2, title('Right target'); end
    
    
    %now for DTTs...
    %For these cases, we will pool the data based on where the obstacle, not which tgt is cued, because presumably they dont know...
%     figure(fn+2);
%     subplot(num_subjects,2,kk+c1); hold on;
     left_obsA_dtt_iang(kk,:) = [iang_dt_left_obsA_bl(kk,:),iang_dt_right_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
%     create_hist_init([iang_dt_left_obsA_bl(kk,:),iang_dt_right_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag) ,purple, 30);
%     if (kk+c1) == 1, title('Leftward obstacle'); end
%     
%     subplot(num_subjects,2,kk+c2); hold on;
     right_obsA_dtt_iang(kk,:) = [iang_dt_right_obsA_bl(kk,:),iang_dt_left_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
%     create_hist_init([iang_dt_right_obsA_bl(kk,:),iang_dt_left_obsA_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag),red_, 30);
%     if (kk+c2) == 2, title('Rightward obstacle'); end
%     
%     %DTT for obstacle B
%     figure(fn+3);
%     subplot(num_subjects,2,kk+c1); hold on;
     left_obsB_dtt_iang(kk,:) = [iang_dt_left_obsB_bl(kk,:),iang_dt_right_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
%     create_hist_init([iang_dt_left_obsB_bl(kk,:),iang_dt_right_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag),purple, 30);
%     if (kk+c1) == 1, title('Leftward obstacle'); end
%     
%     subplot(num_subjects,2,kk+c2); hold on;
     right_obsB_dtt_iang(kk,:) = [iang_dt_right_obsB_bl(kk,:),iang_dt_left_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag);
%     create_hist_init([iang_dt_right_obsB_bl(kk,:),iang_dt_left_obsB_nb(kk,:)] - (bias_dt_null(kk) * dt_bias_flag),red_, 30);
%     if (kk+c2) == 2, title('Rightward obstacle'); end
    
    
%     %make a 2 x2 scatter plot of the average of the median angles for STT
%     %left and right targets for each subject (motor averaging prediction)
%     %vs the observed angle on DTTs
%     figure(fn + 10);
%     %first is obs A Right target
%     subplot(2,2,1); hold on;
%     tmp_ang = median
    
    
    %if kk==2, keyboard; end
    
%     c1 = c1 +1;
%     c2 = c2+1;
    
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
    MT_test = MT_good(exp_seq{3}(end)+1:end); %double check this!
    MT_st_obsA_blk = MT_good( [idx_all.st_left_obsA_bl(kk,:), idx_all.st_right_obsA_bl(kk,:)] );
    MT_st_obsB_blk = MT_good( [idx_all.st_left_obsB_bl(kk,:), idx_all.st_right_obsB_bl(kk,:)] );
    MT_dt_null_test = MT_good( idx_all.dt_null_test(kk,:));
    MT_dt_obsA = MT_good( [idx_all.dt_left_obsA_bl(kk,:), idx_all.dt_right_obsA_bl(kk,:), idx_all.dt_left_obsA_nb(kk,:), idx_all.dt_right_obsA_nb(kk,:)]);
    MT_dt_obsB = MT_good( [idx_all.dt_left_obsB_bl(kk,:), idx_all.dt_right_obsB_bl(kk,:), idx_all.dt_left_obsB_nb(kk,:), idx_all.dt_right_obsB_nb(kk,:)]);
    
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
    MT_st_right_obsA_blk = MT_good( [idx_all.st_right_obsA_bl(kk,:)] );
    MT_st_left_obsA_blk = MT_good( [idx_all.st_left_obsA_bl(kk,:)] );
    
    MT_st_right_obsB_blk = MT_good( [idx_all.st_right_obsB_bl(kk,:)] );
    MT_st_left_obsB_blk = MT_good( [idx_all.st_left_obsB_bl(kk,:)] );
    
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
    
%     figure(fn + 4); 
    
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


%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_test, 'b', []);
%     if c3==1, title('All'); end
%     c3 = c3+1;
% %     xlim([0.4,2]);
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_good(idx_all.st_left_obsA_bl(kk,:)), purple, []); hold on;
%     create_hist_init(MT_good(idx_all.st_right_obsA_bl(kk,:)), red_, []);
%     if c3==2, title('STT obs A (blocked)'); end
%     YL = ylim;
%     mu = nanmean(MT_st_obsA_blk);
%     h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
%     %leg1=legend(h1);
%     %set(leg1,'Location','Best');
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_good(idx_all.st_left_obsB_bl(kk,:)), purple, []); hold on;
%     create_hist_init(MT_good(idx_all.st_right_obsB_bl(kk,:)), red_, []);
%     if c3==3, title('STT obs B (blocked)'); end
%     YL = ylim;
%     mu = nanmean(MT_st_obsB_blk);
%     h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
%     %leg1=legend(h1);
%    %set(leg1,'Location','Best');
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_dt_null_test, 'b', []);
%     p_dtt_null_mt_all(kk) = prctile(MT_dt_null_test,95);
%     if c3==4, title('DTT (null)'); end
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_good([idx_all.dt_left_obsA_bl(kk,:), idx_all.dt_left_obsA_nb(kk,:),]), purple, []); hold on;
%     create_hist_init(MT_good([idx_all.dt_right_obsA_bl(kk,:), idx_all.dt_right_obsA_nb(kk,:)]), red_, []);
%     if c3==5, title('DTT obs A (all)'); end
%     YL = ylim;
%     mu = nanmean(MT_dt_obsA);
%     h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
%     %leg1=legend(h1);
%     %set(leg1,'Location','Best');
%     c3 = c3+1;
%     
%     subplot(num_subjects, 6, c3); hold on;
%     create_hist_init(MT_good([idx_all.dt_left_obsB_bl(kk,:), idx_all.dt_left_obsB_nb(kk,:),]), purple, []); hold on;
%     create_hist_init(MT_good([idx_all.dt_right_obsB_bl(kk,:), idx_all.dt_right_obsB_nb(kk,:)]), red_, []);
%      p_dtt_B_mt_all(kk) = prctile([MT_good([idx_all.dt_left_obsB_bl(kk,:), idx_all.dt_left_obsB_nb(kk,:),]); ...
%          MT_good([idx_all.dt_right_obsB_bl(kk,:), idx_all.dt_right_obsB_nb(kk,:)])],95);
%     if c3==6, title('DTT obs B (all)'); end
%     YL = ylim;
%     mu = nanmean(MT_dt_obsB);
%     h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
%    % leg1=legend(h1);
%     %set(leg1,'Location','Best');
%     c3 = c3+1;
%     
%     
%     figure(333); hold on;
    
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


%     subplot(2,num_subjects,ss); hold on; grid minor;
%     if ~(isempty(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==1)) ))
%         plot(idx_all.st_right_obsA_bl(kk,(opt_indc_st_right_obsA==1)), MT_opt_right_obsA_tmp,'ro');
%     end
%     if ~(isempty(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==1)) ))
%         plot(idx_all.st_left_obsA_bl(kk,(opt_indc_st_left_obsA==1)), MT_opt_left_obsA_tmp,'o','color',purple);
%     end
%     
%    
%     if ~(isempty(idx_all.st_right_obsA_bl((opt_indc_st_right_obsA==0)) ))
%         plot(idx_all.st_right_obsA_bl(kk,(opt_indc_st_right_obsA==0)), MT_ma_right_obsA_tmp,'r.');
%     end
%     if ~(isempty(idx_all.st_left_obsA_bl((opt_indc_st_left_obsA==0)) ))
%         plot(idx_all.st_left_obsA_bl(kk,(opt_indc_st_left_obsA==0)), MT_ma_left_obsA_tmp,'.','color',purple);
%     end
    
    %if ss==1, title('Obs A'); end
    
    
%     subplot(2,num_subjects,ss2); hold on; grid minor;
%     if ~(isempty(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==1)) ))
%         plot(idx_all.st_right_obsB_bl(kk,(opt_indc_st_right_obsB==1)), MT_opt_right_obsB_tmp,'ro');
%     end
%     if ~(isempty(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==1)) ))
%         plot(idx_all.st_left_obsB_bl(kk,(opt_indc_st_left_obsB==1)), MT_opt_left_obsB_tmp,'o','color',purple);
%     end
%     
%     
%     if ~(isempty(idx_all.st_right_obsB_bl((opt_indc_st_right_obsB==0)) ))
%         plot(idx_all.st_right_obsB_bl(kk,(opt_indc_st_right_obsB==0)), MT_ma_right_obsB_tmp,'r.');
%     end
%     if ~(isempty(idx_all.st_left_obsB_bl((opt_indc_st_left_obsB==0)) ))
%         plot(idx_all.st_left_obsB_bl(kk,(opt_indc_st_left_obsB==0)), MT_ma_left_obsB_tmp,'.','color',purple);
%     end

   % if ss==2, title('Obs B'); end
    
%     ss = ss+1;
%     ss2 = ss2+1;
    
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


%% we want to know if there is a significant difference in angle in straight ahead targets between null and obstacle cases

%focus on the blocks that had null and obstacle straight ahead targets in
%the same blocks

%this corresponds to exp_seq 2, 4 and 6
%keyboard;
strt_idx_analysis = [exp_seq{2}, exp_seq{4}, exp_seq{6}]; %blocks we're interested in comparing the cases with

idx_strt_null_sub_all = [];
idx_strt_obs_sub_all = [];


for kk=1:num_subjects
    idx_strt_null_sub_all(kk,:) = strt_idx_analysis( ismember(strt_idx_analysis, idx_all.strt_null(kk,:)));
    
    idx_strt_obs_sub_all(kk,:) = strt_idx_analysis(ismember(strt_idx_analysis, idx_all.strt_obs(kk,:)));       
end

iang_strt_null =nan(num_subjects, size(idx_strt_null_sub_all,2) );
iang_strt_obs =nan(num_subjects, size(idx_strt_obs_sub_all,2) );

for kk2 = 1:num_subjects %can optionally condition to remove/focus only on subjects with obs A data
    iang_strt_null(kk2,:) = initial_angle(kk2, idx_strt_null_sub_all(kk2,:));
    iang_strt_obs(kk2,:) = initial_angle(kk2, idx_strt_obs_sub_all(kk2,:));
    
end

%remove first 5 trials for obs case
num_exclude_strt= 2;
if num_exclude_strt == 0
    iang_strt_obs_tmp = iang_strt_obs;
else 
    iang_strt_obs_tmp = iang_strt_obs(:, num_exclude_strt + 1: end);
end


%reshape it
iang_strt_obs_sub_all = reshape( iang_strt_obs_tmp, numel(iang_strt_obs_tmp), 1);
iang_strt_null_sub_all = reshape( iang_strt_null, numel(iang_strt_null), 1);


% %%show the distributions....
figure; hold on;
create_hist_init_V2(iang_strt_obs_sub_all, ones(3,1)*0.5, 40);
create_hist_init_V2(iang_strt_null_sub_all, [0,0,0], 30);

[h,p] = ttest2(iang_strt_obs_sub_all, iang_strt_null_sub_all);

title(['Straight ahead distribution: obstacle vs no obstacle, p = ', num2str(p)]);



%% Next we need to organize the data according to block
%Eventually, we will sample blocks for each subject in which optimal motions on STTs are made

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
testb = 10; %this is the start of our test block
num_test = length(tpb) - testb +1;


dtt_obsA_right = [];
dtt_obsA_left = [];
dtt_obsB_right = [];
dtt_obsB_left = [];

%%%%%
dt_obsA_right_idx = cell(1,num_subjects);
dt_obsA_left_idx = cell(1,num_subjects);
dt_obsB_right_idx = cell(1,num_subjects);
dt_obsB_left_idx = cell(1,num_subjects);

for qt =1:num_subjects
    dt_obsA_right_idx{qt} = [idx_all.dt_right_obsA_bl(qt,:), idx_all.dt_left_obsA_nb(qt,:)];
    dt_obsA_left_idx{qt} = [idx_all.dt_left_obsA_bl(qt,:), idx_all.dt_right_obsA_nb(qt,:)];
    
    dt_obsB_right_idx{qt} = [idx_all.dt_right_obsB_bl(qt,:), idx_all.dt_left_obsB_nb(qt,:)];
    dt_obsB_left_idx{qt} = [idx_all.dt_left_obsB_bl(qt,:), idx_all.dt_right_obsB_nb(qt,:)];
end
% dt_obsA_right_idx = [idx_all.dt_right_obsA_bl, idx_all.dt_left_obsA_nb];
% dt_obsA_left_idx = [idx_all.dt_left_obsA_bl, idx_all.dt_right_obsA_nb];

% dt_obsB_right_idx = [idx_all.dt_right_obsB_bl, idx_all.dt_left_obsB_nb];
% dt_obsB_left_idx = [idx_all.dt_left_obsB_bl, idx_all.dt_right_obsB_nb];

%indicator for which blocks to sample from which subjects
opt_block_sample = nan(num_subjects,num_test);


block_st_right_obsA_bl_idx = cell(num_subjects,num_test);
block_st_left_obsA_bl_idx = cell(num_subjects,num_test);

block_st_right_obsB_bl_idx = cell(num_subjects,num_test);
block_st_left_obsB_bl_idx = cell(num_subjects,num_test);

block_dt_obsA_right_idx = cell(num_subjects,num_test);
block_dt_obsA_left_idx = cell(num_subjects,num_test);
block_dt_obsB_right_idx = cell(num_subjects,num_test);
block_dt_obsB_left_idx = cell(num_subjects,num_test);


for jj=1:num_test
    for kk=1:num_subjects
        
        current_block = jj + (testb) - 1;
        current_idx = block_ind{current_block};
        
        %get the indices we want
        %B needs to be the current_idx
        
        block_st_right_obsA_bl_idx{kk,jj} = idx_all.st_right_obsA_bl(kk,ismember(idx_all.st_right_obsA_bl(kk,:), current_idx));
        block_st_left_obsA_bl_idx{kk,jj} = idx_all.st_left_obsA_bl(kk,ismember(idx_all.st_left_obsA_bl(kk,:), current_idx));
        
        block_st_right_obsB_bl_idx{kk,jj} = idx_all.st_right_obsB_bl(kk,ismember(idx_all.st_right_obsB_bl(kk,:), current_idx));
        block_st_left_obsB_bl_idx{kk,jj} = idx_all.st_left_obsB_bl(kk,ismember(idx_all.st_left_obsB_bl(kk,:), current_idx));
        
        block_dt_obsA_right_idx{kk,jj} = dt_obsA_right_idx{kk}(ismember( dt_obsA_right_idx{kk}, current_idx));
        block_dt_obsA_left_idx{kk,jj} = dt_obsA_left_idx{kk}(ismember( dt_obsA_left_idx{kk}, current_idx));
        
        block_dt_obsB_right_idx{kk,jj} = dt_obsB_right_idx{kk}(ismember( dt_obsB_right_idx{kk}, current_idx));
        block_dt_obsB_left_idx{kk,jj} = dt_obsB_left_idx{kk}(ismember( dt_obsB_left_idx{kk}, current_idx));

        
        %figure out if during certain blocks, a subject produced at least 70% optimal movements on STTs
        opt_idx = optimal_indicator_all(kk,:);
        
        block_opt_idx_obsA_right = opt_idx(block_st_right_obsA_bl_idx{kk,jj});
        block_opt_idx_obsA_left = opt_idx(block_st_left_obsA_bl_idx{kk,jj});
        
        block_opt_idx_obsB_right = opt_idx(block_st_right_obsB_bl_idx{kk,jj});
        block_opt_idx_obsB_left = opt_idx(block_st_left_obsB_bl_idx{kk,jj});
        
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

%% figure(555); hold on;
%now lets get the data one want from the good blocks

% for p = 1:num_subjects%6
%     for q=1:num_test
%         
%         %sample the blocks we want
%         if opt_block_sample(p,q)
%             %get the angle!!! (and also trajectories???)
%             
%             for i=1:length(block_dt_obsA_right_idx{q})
%                 dtt_obsA_right = [dtt_obsA_right, initial_angle(p,block_dt_obsA_right_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
%                 
%                  %%%check some random trials
% %                  if initial_angle(p,block_dt_obsA_right_idx{q}(i)) <= -60
% %                      keyboard;
% %                  end
%                 
%             end
%             
%             for i=1:length(block_dt_obsA_left_idx{q})
%                 dtt_obsA_left = [dtt_obsA_left, initial_angle(p,block_dt_obsA_left_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
%             end
%             
%             for i=1:length(block_dt_obsB_right_idx{q})
%                 dtt_obsB_right = [dtt_obsB_right, initial_angle(p,block_dt_obsB_right_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
%             end
%             
%             for i=1:length(block_dt_obsB_left_idx{q})
%                 dtt_obsB_left = [dtt_obsB_left, initial_angle(p,block_dt_obsB_left_idx{q}(i)) - (bias_dt_null(p) * dt_bias_flag) ];
%             end
%             
% %             if p==6,
% %                 figure(555); hold on;
% %                 for zz=1:length(block_dt_obsA_left_idx{q})
% %                     plot(initial_mov{p,block_dt_obsA_left_idx{q}(zz)}(:,1), initial_mov{p,block_dt_obsA_left_idx{q}(zz)}(:,2) );
% %                 end
% %             end
%             
%             %if p==6 & q==3, keyboard; end
%             
%             
%         end
%     end
% end
% 
% 
% figure(fn+ 6);
% subplot(2,2,1); hold on;
% create_hist_init(dtt_obsA_right, grey, []);
% YL = ylim;
% mu = nanmean(dtt_obsA_right);
% h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% %leg1=legend(h1);
% %set(leg1,'Location','Best');
% title('Right obstacle');
% ylabel('obs A');
% subplot(2,2,2); hold on;
% create_hist_init(dtt_obsA_left, grey, []);
% title('Left Obstacle');
% YL = ylim;
% mu = nanmean(dtt_obsA_left);
% h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% % leg1=legend(h1);
% % set(leg1,'Location','Best');
% 
% 
% subplot(2,2,3); hold on;
% ylabel('obs B');
% create_hist_init(dtt_obsB_right, grey, []);
% YL = ylim;
% mu = nanmean(dtt_obsB_right);
% h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% % leg1=legend(h1);
% % set(leg1,'Location','Best');
% 
% subplot(2,2,4); hold on;
% create_hist_init(dtt_obsB_left, grey, []);
% YL = ylim;
% mu = nanmean(dtt_obsB_left);
% h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% % leg1=legend(h1);
% % set(leg1,'Location','Best');



%% optionally filter the data based on iqr!!!!

%NOTE: We might want to do this earlier on for the DTT null (baseline) data as well
iqr_thresh = 3; %filter data based on 3 iqrs away

%idx_left_obsA_dtt = []

if iqr_filter_flag
    
    %filter obs DTT data
    [right_obsA_dtt_iang, right_obsA_dtt_iqr_idx] = filter_iqr_0415a(right_obsA_dtt_iang, iqr_thresh);
    [left_obsA_dtt_iang, left_obsA_dtt_iqr_idx] = filter_iqr_0415a(left_obsA_dtt_iang, iqr_thresh);
    
    %repeat for baseline DTT data
    [iang_dt_bln, iang_dt_bln_iqr_idx] = filter_iqr_0415a(iang_dt_bln, iqr_thresh);
    
    %lets see if the iqr captures the inside movements on stt trials with obstacles
    [iang_st_left_obsA_bl, iang_st_left_obsA_bl_iqr_idx] = filter_iqr_0415a(iang_st_left_obsA_bl, iqr_thresh);
    [iang_st_right_obsA_bl, iang_st_right_obsA_bl_iqr_idx] = filter_iqr_0415a(iang_st_right_obsA_bl, iqr_thresh);
    
    %try no block case as well, even though it probably doesnt matter much
    [iang_st_left_obsA_nb, iang_st_left_obsA_nb_iqr_idx] = filter_iqr_0415a(iang_st_left_obsA_nb, iqr_thresh);
    [iang_st_right_obsA_nb, iang_st_right_obsA_nb_iqr_idx] = filter_iqr_0415a(iang_st_right_obsA_nb, iqr_thresh);
    
    %repeat for baseline STT data
    [iang_st_left_bln, iang_st_left_iqr_idx] = filter_iqr_0415a(iang_st_left_bln, iqr_thresh);
    [iang_st_right_bln, iang_st_right_iqr_idx] = filter_iqr_0415a(iang_st_right_bln, iqr_thresh);
    [iang_st_strt_bln, iang_st_strt_iqr_idx] = filter_iqr_0415a(iang_st_strt_bln, iqr_thresh);
    
    %count what we filtered
    num_filt_iqr_dt_trials = sum(sum(left_obsA_dtt_iqr_idx,1)) + sum(sum(right_obsA_dtt_iqr_idx,1));
    frac_filt_iqr_dt_trials = num_filt_iqr_dt_trials / ( numel(left_obsA_dtt_iqr_idx) + numel(right_obsA_dtt_iqr_idx) );
    
    num_filt_iqr_st_trials = sum(sum(iang_st_left_obsA_bl_iqr_idx,1)) + sum(sum(iang_st_right_obsA_bl_iqr_idx,1)) +...
        sum(sum(iang_st_left_obsA_nb_iqr_idx,1)) + sum(sum(iang_st_right_obsA_nb_iqr_idx,1)) +...
        sum(sum(iang_dt_bln_iqr_idx,1)) + sum(sum(iang_st_left_iqr_idx,1)) + sum(sum(iang_st_right_iqr_idx,1)) +...
        sum(sum(iang_st_strt_iqr_idx,1));
    frac_filt_iqr_st_trials = num_filt_iqr_st_trials / ( numel(iang_st_left_obsA_bl_iqr_idx) + numel(iang_st_right_obsA_bl_iqr_idx) + ...
        numel(iang_st_left_obsA_nb_iqr_idx) + numel(iang_st_right_obsA_nb_iqr_idx) + numel(iang_dt_bln_iqr_idx) +...
        numel(iang_st_left_iqr_idx) + numel(iang_st_right_iqr_idx) + numel(iang_st_left_iqr_idx) );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %check to see if iqr caught the bad st trials
    
    %     tmp_data1 =iang_st_right_obsA_bl_iqr_idx;
    %     tmp_data2 = inside_indc.obsA_R;
    %
    %     cc1 = reshape(tmp_data1, [size(tmp_data1,1)*size(tmp_data1,2), 1]);
    %     cc2 = reshape(tmp_data2, [size(tmp_data2,1)*size(tmp_data2,2), 1]);
    %
    %     find(cc1)
    %     find(cc2)
    
    
    %get the indices of bad trials corresponding to the overall experiment so
    %we can filter these out in other scripts
    sub_iqr_bad_trial_idx = cell(num_subjects, 1);
    
    for np=1:num_subjects
        
        %get the right indices for the dt trials
        idx_left_obsA_sub_tmp = [idx_all.dt_left_obsA_bl(np,:), idx_all.dt_right_obsA_nb(np,:)];
        idx_right_obsA_sub_tmp = [idx_all.dt_right_obsA_bl(np,:), idx_all.dt_left_obsA_nb(np,:)];
        
        sub_iqr_bad_trial_idx{np} = [idx_all.st_left_obsA_bl(np, find(iang_st_left_obsA_bl_iqr_idx(np,:))), ...
            idx_all.st_right_obsA_bl(np, find(iang_st_right_obsA_bl_iqr_idx(np,:))), idx_all.st_left_obsA_nb(np, find(iang_st_left_obsA_nb_iqr_idx(np,:))),...
            idx_all.st_right_obsA_nb(np, find(iang_st_right_obsA_nb_iqr_idx(np,:))), idx_left_obsA_sub_tmp(find(left_obsA_dtt_iqr_idx(np,:))),...
            idx_all.st_left_bln(np, find(iang_st_left_iqr_idx(np,:))), idx_all.st_right_bln(np, find(iang_st_right_iqr_idx(np,:))),...
            idx_all.st_strt_bln(np, find(iang_st_strt_iqr_idx(np,:))), idx_all.dt_bln(np, find(iang_dt_bln_iqr_idx(np,:))),...
            idx_right_obsA_sub_tmp(find(right_obsA_dtt_iqr_idx(np,:)))];
    end
    
    save('sub_iqr_bad_trial_idx.mat','sub_iqr_bad_trial_idx');
    
end

%[iang_dt_right_obsA_bl(kk,:),iang_dt_left_obsA_nb(kk,:)]
% [iang_dt_right_obsA_bl(kk,:),iang_dt_left_obsA_nb(kk,:)]
% idx_all.dt_left_obsA_bl

%% MOVEMENT ANGLE HISTOGRAMS

%     iang_dt_bln(qs,:) = initial_angle(qs, idx_all.dt_bln(qs,:));
%     iang_st_left_bln(qs,:) = initial_angle(qs,idx_all.st_left_bln(qs,:));
%     iang_st_right_bln(qs,:) = initial_angle(qs,idx_all.st_right_bln(qs,:));

right_obsA_dtt_iang_all = reshape(right_obsA_dtt_iang,numel(right_obsA_dtt_iang),1);
left_obsA_dtt_iang_all = reshape(left_obsA_dtt_iang,numel(left_obsA_dtt_iang),1) ;

null_dtt_iang_all = reshape(iang_dt_bln, numel(iang_dt_bln), 1);

figure; hold on; title('distribution of 2-target angles across all subjects');
create_hist_init(right_obsA_dtt_iang_all, red_, 40);
create_hist_init(left_obsA_dtt_iang_all, light_blue, 40);
%create_hist_init_sc(null_dtt_iang_all, grey, 40);

YL = ylim;
mu1 = nanmean(right_obsA_dtt_iang_all);
mu2 = nanmean(left_obsA_dtt_iang_all);
mu3 = nanmean(null_dtt_iang_all);

h1=plot([mu1, mu1], [YL(1), YL(2)], 'color', red_, 'linestyle','--','linewidth',1.5,'Displayname',num2str(mu1));
h2=plot([mu2, mu2], [YL(1), YL(2)], 'color', light_blue, 'linestyle','--','linewidth',1.5,'Displayname',num2str(mu2));

%h3=plot([mu3, mu3], [YL(1), YL(2)], 'color', grey, 'linestyle',':','linewidth',1.5,'Displayname',num2str(mu3));
% leg1=legend(h1);
% set(leg1,'Location','Best');
 xlim([-40,40]);

%do the same for the STT data
%right and left ->both with and without obstacles
%NOTE: this doesn't include trials in which the obstacle was present but
%didnt directly block movement to the target

figure; hold on; title('distribution of 1-target angles across all subjects');


st_left_bln_all = reshape(iang_st_left_bln, numel(iang_st_left_bln), 1);
st_right_bln_all = reshape(iang_st_right_bln, numel(iang_st_right_bln), 1);

st_left_obsA_bl_all = reshape(iang_st_left_obsA_bl, numel(iang_st_left_obsA_bl),1);
st_right_obsA_bl_all = reshape(iang_st_right_obsA_bl, numel(iang_st_right_obsA_bl),1);

create_hist_init(st_left_obsA_bl_all, light_blue, 40);
create_hist_init(st_right_obsA_bl_all, red_, 40);

% create_hist_init_sc(st_left_bln_all, light_blue, 40);
% create_hist_init_sc(st_right_bln_all, red_, 40);

%xlim([-60,60]);

YL = ylim;

mu1 = nanmean(st_left_obsA_bl_all);
mu2 = nanmean(st_right_obsA_bl_all);
% mu3 = nanmean(st_right_bln_all);
% mu4 = nanmean(st_left_bln_all);

h1=plot([mu1, mu1], [YL(1), YL(2)], 'color', light_blue, 'linestyle','--','linewidth',1.5);
h2=plot([mu2, mu2], [YL(1), YL(2)], 'color', red_, 'linestyle','--','linewidth',1.5);

% h3=plot([mu3, mu3], [YL(1), YL(2)], 'color', red_, 'linestyle',':','linewidth',1.5);
% h4=plot([mu4, mu4], [YL(1), YL(2)], 'color', light_blue, 'linestyle',':','linewidth',1.5);

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

%%


avg_pred_obsA_right = ( (nanmedian(iang_st_left_obsA_nb,2) - bias_st_left *st_bias_flag) + ...
    (nanmedian(iang_st_right_obsA_bl,2) -bias_st_right *st_bias_flag) ) /2;
%avg_obsA_right =  (nanmedian(iang_st_right_obsA_bl,2) -bias_st_right *st_bias_flag)+30;
dtt_obsA_right_all = nanmedian(right_obsA_dtt_iang,2);

avg_pred_obsA_left = ( (nanmedian(iang_st_right_obsA_nb,2) - bias_st_right *st_bias_flag) + ...
    (nanmedian(iang_st_left_obsA_bl,2) - bias_st_left *st_bias_flag) ) /2;
%avg_obsA_left = nanmedian(iang_st_left_obsA_bl,2) - bias_st_left *st_bias_flag - 30;
dtt_obsA_left_all = nanmedian(left_obsA_dtt_iang,2);


figure; hold on;
xlabel('MA prediction');
ylabel('Observed angle on DTT');

%should we combine left and right data?
avg_pred_obsA_LR = [-1*avg_pred_obsA_right; avg_pred_obsA_left];
dtt_combined_all = [-1*dtt_obsA_right_all; dtt_obsA_left_all]; 

%get the correlation coeffcient
[r_ma_pred, p_ma_pred] = corr(avg_pred_obsA_LR, dtt_combined_all);
[r_ma_pred_obsA_right, p_ma_pred_obsA_right] = corr(-1*avg_pred_obsA_right, -1*dtt_obsA_right_all);
[r_ma_pred_obsA_left, p_ma_pred_obsA_left] = corr(avg_pred_obsA_left, dtt_obsA_left_all);


title(['correlation is r = ', num2str(r_ma_pred), ',  p = ', num2str(p_ma_pred)]);

% plot(avg_pred_obsA_LR, dtt_combined_all, '.', 'markersize', 15);
plot(avg_pred_obsA_right, dtt_obsA_right_all, 'r.', 'markersize', 15);
plot(avg_pred_obsA_left, dtt_obsA_left_all, 'b.', 'markersize', 15);

%add the unity line to show the MA prediction
plot([-30:30], [-30:30], 'k--', 'linewidth', 1);

ylim([-30,30]);

%add standard error maybe???
%----> cant, we dont have equal number of trials (pairs)!


%in this case each trial is subtracted by the bias
% avg_pred_obsA_right_sub = ( iang_st_left_obsA_nb - ((bias_st_left * st_bias_flag) * ones(1,size(iang_st_left_obsA_nb,2))) )  + ... 
%     ( iang_st_right_obsA_bl - ((bias_st_right * st_bias_flag) * ones(1,size(iang_st_right_obsA_bl,2))) ) /2;


%%

%phiA = atand(2/9); %how much the obstacle sticks out, and where it is frist seen along its axis
phiA = atand(2/(obs_locy - sind(60)));

% %cut it here
% left_obsA_dtt_iang_cut = left_obsA_dtt_iang(: ,num_cutoff+1:end);
% right_obsA_dtt_iang_cut = right_obsA_dtt_iang(: ,num_cutoff+1:end);

%look at the last x # of trials to see if they settle around an angle
num_end = 0; %look at the last x number of trials

% s_AL = st_block_dat.obsA_L;
% s_AR = st_block_dat.obsA_R;

if num_end >0
    s_AL = st_block_dat.obsA_L(:, end - num_end+1:end);
    s_AR = st_block_dat.obsA_R(:, end - num_end+1:end);
    
    left_obsA_dtt_iang_cut = left_obsA_dtt_iang(: ,end - num_end+1:end);
    right_obsA_dtt_iang_cut = right_obsA_dtt_iang(: ,end - num_end+1:end);
    
else
    
    s_AL = st_block_dat.obsA_L;
    s_AR = st_block_dat.obsA_R;
    
    left_obsA_dtt_iang_cut = left_obsA_dtt_iang;
    right_obsA_dtt_iang_cut = right_obsA_dtt_iang;
end
        
%%

%Lets plot safety margin on stt vs observed angle on dtt, but focus on magnitude here
%can also consider subtracting dtt data by the amount the obstacle sticks out

% figure; hold on;
% xlabel('Safety margin on STT');
% ylabel('Observed angle on DTT');
% % 
 dtt_iang_comb_mag = [nanmedian( left_obsA_dtt_iang*-1,2); nanmedian(right_obsA_dtt_iang,2)];
% % % stt_SM_comb_mag = [nanmean(st_block_dat.obsA_L,2); nanmean( st_block_dat.obsA_R * -1,2)]; %safety margin
% % 
% dt_leftA_sub = nanmean(left_obsA_dtt_iang_cut*-1,2);
% st_leftA_SM_sub = nanmean(s_AL,2);
% 
% dt_rightA_sub = nanmean(right_obsA_dtt_iang_cut,2);
% st_rightA_SM_sub = nanmean(s_AR*-1,2);
% 
% 
% plot(st_leftA_SM_sub,dt_leftA_sub ,'.');
% plot(st_rightA_SM_sub,dt_rightA_sub ,'.');

%%%corr(st_leftA_SM_sub, dt_leftA_sub)
%%%corr(st_rightA_SM_sub, dt_rightA_sub)


%% Lets see if baseline variability on DTT can predict size of safety margins....

%estimate baseline variability for each subject

baseline_var_dt_bln = nan(num_subjects, 1);
baseline_var_st_right = nan(num_subjects, 1);
baseline_var_st_left = nan(num_subjects, 1);
baseline_var_st_strt = nan(num_subjects, 1);

test_var_dt_right = nan(num_subjects, 1);
test_var_dt_left = nan(num_subjects, 1);
test_var_dt_comb = nan(num_subjects, 1);

pvar_all = nan(num_subjects,1);
pse_all = nan(num_subjects, 1);

iqr_dt_test_sub = nan(2, num_subjects);
std_dt_test_sub = nan(2, num_subjects);

%we will also want the std of obstacle-free and obstacle-present 1-target
%trials (use only test block data)
std_st_test_obs_free = nan(2,num_subjects); %left and right - could combine straight as well?
std_st_test_obs_present = nan(2,num_subjects);

% st_left_bln: [26×20 double]
% st_right_bln: [26×20 double]
% st_strt_bln: [26×20 double]
% right_obsA_dtt_iang

for kq2 = 1:num_subjects
    
    iang_st_right = initial_angle(kq2, idx_all.st_right_bln(kq2,:)) - (bias_st_right(kq2) * st_bias_flag);
    iang_st_left = initial_angle(kq2, idx_all.st_left_bln(kq2,:)) - (bias_st_left(kq2) * st_bias_flag);
    iang_st_strt = initial_angle(kq2, idx_all.st_strt_bln(kq2,:)) - (bias_st_strt(kq2) * st_bias_flag);
    
    %calculate the variance
    baseline_var_st_right(kq2) = nanvar(iang_st_right,0,2);
    baseline_var_st_left(kq2) = nanvar(iang_st_left,0,2);
    baseline_var_st_strt(kq2) = nanvar(iang_st_strt,0,2);
    baseline_var_dt_bln(kq2) = nanvar(iang_dt_bln(kq2, :),0,2);
    
    test_var_dt_right(kq2) = nanvar(right_obsA_dtt_iang_cut(kq2,:),0,2);
    test_var_dt_left(kq2) = nanvar(left_obsA_dtt_iang_cut(kq2,:),0,2);
    test_var_dt_comb(kq2) = nanvar([right_obsA_dtt_iang_cut(kq2,:), left_obsA_dtt_iang_cut(kq2,:)],0,2);
    
    iqr_dt_test_sub(1,kq2) = iqr(left_obsA_dtt_iang_cut(kq2,:));
    iqr_dt_test_sub(2,kq2) = iqr(right_obsA_dtt_iang_cut(kq2,:));
    
    %test during 1-target trials
    std_st_test_obs_free(1,kq2) = nanstd(iang_st_left_obsA_nb(kq2,:));
    std_st_test_obs_free(2,kq2) = nanstd(iang_st_right_obsA_nb(kq2,:));
    std_st_test_obs_present(1,kq2) = nanstd(s_AL(kq2,:));
    std_st_test_obs_present(2,kq2) = nanstd(s_AR(kq2,:));
    
    %collect individual subject pvalues to make sure the variance isnt
    %different on an individual subject level
    %if kq2==3, keyboard; end
    [pvar_all(kq2), pse_all(kq2)] = permtest_var_1210_2019a(right_obsA_dtt_iang_cut(kq2,:), [left_obsA_dtt_iang_cut(kq2,:),1,2,3,4,5], 1, 1e3);
    %[~,pvar_all(kq2)] = vartest2(right_obsA_dtt_iang(kq2,:), left_obsA_dtt_iang(kq2,:));  
end

%calculate the std as well
std_dt_test_sub(1,:) = sqrt(test_var_dt_left);
std_dt_test_sub(2,:) = sqrt(test_var_dt_right);
std_dt_test_sub_comb = (sqrt(test_var_dt_left) + sqrt(test_var_dt_right))/2;

%smallest p value can only be the resolution of number of samples (default is 1e3)
pvar_all(pvar_all==0) = 1e-3;

figure; hold on; 
xlabel('Baseline std');
ylabel('Observed angle on DTT');

period_flag = 0;
early_idx_tmp = 15; %15 trials is the first quarter

if period_flag
    qtmp_left_obsA_dtt_iang = left_obsA_dtt_iang(:,1:early_idx_tmp);
    qtmp_right_obsA_dtt_iang = right_obsA_dtt_iang(:,1:early_idx_tmp);
else
    qtmp_left_obsA_dtt_iang = left_obsA_dtt_iang;
    qtmp_right_obsA_dtt_iang = right_obsA_dtt_iang;
end

plot([baseline_var_dt_bln; baseline_var_dt_bln], dtt_iang_comb_mag, '.');
plot([baseline_var_dt_bln], nanmedian( qtmp_left_obsA_dtt_iang*-1,2), 'b.');
plot([baseline_var_dt_bln], nanmedian(qtmp_right_obsA_dtt_iang,2), 'r.');

corr([baseline_var_dt_bln; baseline_var_dt_bln], dtt_iang_comb_mag)
corr([baseline_var_dt_bln], nanmedian( qtmp_left_obsA_dtt_iang*-1,2))
corr([baseline_var_dt_bln], nanmedian( qtmp_right_obsA_dtt_iang,2))

%%%calculate the ratio of 2-target trial angle std to 1-target trial angle
%%%std averaged over all three target directions
dt_bln_std_avg = nanmedian(sqrt(baseline_var_dt_bln));

st_bln_std_right_avg = nanmedian(sqrt(baseline_var_st_right));
st_bln_std_left_avg = nanmedian(sqrt(baseline_var_st_left));
st_bln_std_strt_avg = nanmedian(sqrt(baseline_var_st_strt));

st_bln_std_avg = (st_bln_std_right_avg + st_bln_std_left_avg + st_bln_std_strt_avg) / 3;

%std_ratio = dt_bln_std_avg/st_bln_std_avg;

std_st_all = median([sqrt(baseline_var_st_right), sqrt(baseline_var_st_left), sqrt(baseline_var_st_strt)],2);
std_dt_all = sqrt(baseline_var_dt_bln);

%load in the experiment 2a data
exp2a = load('exp2a_corr_data.mat');

%estimate the slope by combining data from 2a and 2b
std_ratio_comb = [std_ratio_all; exp2a.std_ratio_all];

final_std_ratio = median(std_ratio_comb); %2.1459

%linspace(5,95,5)

%% plot histogram of p values for left vs right variance
figure; hold on;
histogram(pvar_all,24,'facealpha',0.75, 'edgecolor','k');
ylim([0, length(pvar_all(pvar_all<=0.05))+0]);
vline(0.05, 'k--');

title('comparison of left vs right 2TT variance');
xlabel('p-value');
ylabel('# of participants');

%get the subjects that supposedly have a difference and make sure we're not dealing with outliers
%we will show subjects that arent just less than 0.05, but are also borderline based on an estimate 
%of the resampling noise from the permutation test -- we can just use the upper 95% CI
hd_sub_id = find(pvar_all<(0.05 + 1.96/2*median(pse_all(pse_all~=0)))); %"high difference" idk what else to call the damn thing

%save number of "bad" samples for these subjects
nb_hd_sub = nan(2, length(hd_sub_id));

%remove_sample flag determines if we should remove samples based on iqr for
%these 5 subjects
remove_sample_flag = 1;
[~,num_dtt_samples] = size(right_obsA_dtt_iang);
right_dtt_filt = nan(length(hd_sub_id), num_dtt_samples);
left_dtt_filt = nan(length(hd_sub_id), num_dtt_samples);

figure;
for k = 1:length(hd_sub_id)
    subplot(ceil(length(hd_sub_id)/3), 3, k); hold on;
    %optionally show the histogram with mean subtracted to focus on the spread of the data
%     histogram(right_obsA_dtt_iang(hd_sub_id(k),:) - nanmean(right_obsA_dtt_iang(hd_sub_id(k),:)),20,'facecolor', 'r','facealpha',0.5, 'edgecolor','k');
%     histogram(left_obsA_dtt_iang(hd_sub_id(k),:) - nanmean(left_obsA_dtt_iang(hd_sub_id(k),:)),20,'facecolor', 'b','facealpha',0.5, 'edgecolor','k');
    histogram(right_obsA_dtt_iang_cut(hd_sub_id(k),:),20,'facecolor', 'r','facealpha',0.5, 'edgecolor','k');
    histogram(left_obsA_dtt_iang_cut(hd_sub_id(k),:),20,'facecolor', 'b','facealpha',0.5, 'edgecolor','k');
    
    xlabel('angle');
    ylabel('count');    
    
    %for these subjects, count the number of sample above/below x times the iqr
    [left_dtt_filt(k,:), b1_tmp] = filter_sub_iqr_0527_2020a(left_obsA_dtt_iang_cut(hd_sub_id(k),:), 2);
    [right_dtt_filt(k,:), b2_tmp] = filter_sub_iqr_0527_2020a(right_obsA_dtt_iang_cut(hd_sub_id(k),:), 2);
    
    nb_hd_sub(1,k) = sum(b1_tmp);
    nb_hd_sub(2,k) = sum(b2_tmp);
    
    %now redo comparison of p values
    [p_tmp, ~] = permtest_var_1210_2019a(right_dtt_filt(k,:), left_dtt_filt(k,:), 0);

%     title(['p = ', num2str(pvar_all(hd_sub_id(k))), ' w/out filter & p = ', num2str(p_tmp), ' w/ filter, n = ',...
%         num2str(sum(nb_hd_sub(:,k))), ' samples removed']);
%     title(['p = ', num2str(pvar_all(hd_sub_id(k))), ' w/out filter & p = ', num2str(p_tmp), ' w/ filter']);
    %title(['p = ', num2str(pvar_all(hd_sub_id(k)))]);
    title({['right mean = ', num2str(nanmean(right_obsA_dtt_iang_cut(hd_sub_id(k),:))) ' left mean = ' num2str(nanmean(left_obsA_dtt_iang_cut(hd_sub_id(k),:)))],...
        ['right median = ', num2str(nanmedian(right_obsA_dtt_iang_cut(hd_sub_id(k),:))) ' left median = ' num2str(nanmedian(left_obsA_dtt_iang_cut(hd_sub_id(k),:)))],...
        ['right std = ', num2str(std_dt_test_sub(2,hd_sub_id(k))), ' left std = ' num2str(std_dt_test_sub(1,hd_sub_id(k)))],...
        ['right se = ', num2str(std_dt_test_sub(2,hd_sub_id(k))/sqrt(sum(~isnan(right_obsA_dtt_iang_cut(hd_sub_id(k),:))))),...
        ' left se = ' num2str(std_dt_test_sub(1,hd_sub_id(k))/sqrt(sum(~isnan(left_obsA_dtt_iang_cut(hd_sub_id(k),:)))))]});
end
%% lets look at the instantaneous rwd rate on blocked st trials as well as dt trials

st_block_test_rwd_all = struct(obs_fields{1}, [], obs_fields{2}, [] );
st_block_test_rwd_all.comb = [];

%use sub_rwd_all
%dont bother zero-ing when test phase starts, since theres no reason to
%expect that changes throughout the experiment shouldnt affect certain
%trials more, other than blocking the target directly

rwd_agg = cumsum(sub_rwd_all,1)./repmat([1:size(sub_rwd_all,1)],num_subjects,1)';

%%get the indices correspinding to when participants are in a test phase,
%%but only for obs A

for k=1:num_subjects
    idx_tmp = find(info_all.is_obsA_all(k,:));
    testA_idx(k,:) = idx_tmp(idx_tmp>exp_seq{6}(end));
    
    idx_st_right_testA_tmp = find(info_all.is_obsA_all(k,:) .* info_all.isRight_all(k,:) .* info_all.ST_idx_all(k,:));
    idx_st_left_testA_tmp = find(info_all.is_obsA_all(k,:) .* info_all.isLeft_all(k,:) .* info_all.ST_idx_all(k,:));

    testA_st_right_idx(k,:) = find(idx_st_right_testA_tmp>exp_seq{6}(end));
    testA_st_left_idx(k,:) = find(idx_st_left_testA_tmp>exp_seq{6}(end));
    
    rwd_testA(:,k) = sub_rwd_all(testA_idx(k,:),k);
end

rwd_testA_agg = cumsum(rwd_testA,1)./repmat([1:size(rwd_testA,1)],num_subjects,1)';


% for qz =1:num_subjects
%     for w1=1:length(obs_fields)
%        rwd_trials_tmp = rwd_label(qz, stt_block_idx.(obs_fields{w1}) );
%        st_block_test_rwd_all.(obs_fields{w1})(qz) = sum(rwd_trials_tmp) / length(rwd_trials_tmp) * 100;
%     end
% end
% 
% %look at blocked trials in aggregate
% for qz2=1:num_subjects
%     st_block_test_rwd_all.comb(qz2) = sum(rwd_label(qz2, st_block_idx_comb)) / length(st_block_idx_comb) *100;
% end

%plot the scores per block
% figure; hold on;
% for qz5 = 1:num_subjects
%    plot( info_all.Score_sub_all(qz5,:));
% end
% 
% xlabel('exp block');
% ylabel('Score');

%could probably also make a time series of the feedback given (check it
%early on especially...)

%% do MA and PO predictions

%close all;
% figure(fn + 10); hold on;
% title('left and right data, angle = obstacle axis, no bias subtraction');
% 
% %start with obs A, left
% display_PO_0519_2020a( s_AL, left_obsA_dtt_iang_cut, -1, phiA, [], std_ratio); %need input for population averaged ratio
% 
% %obs A, right
% display_PO_0519_2020a( s_AR, right_obsA_dtt_iang_cut, 1, phiA, [], std_ratio);

%%
figure; hold on;
title('exp 2b predictions');

%subtracting with NaNs will cause us to lose samples! Take mean/median first
% sm_diff = (nanmean(s_AL,2) - nanmean(s_AR,2))/2;
% dt_diff = (nanmean(left_obsA_dtt_iang_cut,2) - nanmean(right_obsA_dtt_iang_cut,2))/2;

sm_diff = (-nanmean(s_AL,2) + nanmean(s_AR,2))/2;
dt_diff = (-nanmedian(left_obsA_dtt_iang_cut,2) + nanmedian(right_obsA_dtt_iang_cut,2))/2;
%dt_diff = nanmedian(-left_obsA_dtt_iang_cut - nanmedian(right_obsA_dtt_iang_cut,2))/2;

stmp = (s_AR - s_AL)/2;
stmp = bsxfun(@minus, stmp, nanmedian(stmp));
std_st_test_obs_present_comb = nanstd(stmp, 0, 2);

dtmp = (right_obsA_dtt_iang_cut - left_obsA_dtt_iang_cut)/2;
dtmp = bsxfun(@minus, dtmp, nanmedian(dtmp));
std_dt_test_obs_present_comb = nanstd(dtmp, 0, 2);

%calculate the population averaged ratio by first calculating the mean SDs
%and then calculting the ratios. This avoids bias induced from dividing by
%a normally distributed variable (mean will always be greater than 1)
obs_ratio_gm = mean(std_dt_test_sub_comb)/mean(mean(std_st_test_obs_present,1)); %about 1.003

%calculate each subject's variability
obs_ratio_sub = std_dt_test_sub_comb./mean(std_st_test_obs_present,1)';

s2_tmp1 = std_dt_test_obs_present_comb;
s2_tmp2 = std_dt_test_sub_comb;

s1_tmp1 = std_st_test_obs_present_comb;
s1_tmp2 = mean(std_st_test_obs_present,1)';

%we will want to factor the bias from STT non-obstructed movements in the MA predictions
obs_nb_bias.left = (nanmedian(iang_st_left_obsA_nb,2)) - 30;
obs_nb_bias.right =(nanmedian(iang_st_right_obsA_nb,2)) + 30;
obs_nb_bias.comb = (obs_nb_bias.left + obs_nb_bias.right)/2;

[~, pmodel, smodel, r2_model, corr_pred, p_corr_pred, dt_diff_ma_shifted, dt_diff_po_shifted, po_pred, ma_pred] = ...
    display_PO_diff_flipped_0623_2020a( sm_diff, dt_diff, 1, phiA, [], 0.66, 0.71, obs_nb_bias, obs_ratio_sub, obs_ratio_gm); %refined slopes are from Exp 2a

%% estimate what the slope would be in exp2b data (but dont use it for the predictions)
[b_pref_ma_check, b_pref_po_check] = determine_model_weights(sm_diff, dt_diff_ma_shifted, dt_diff_po_shifted, [], [], obs_ratio_sub); %avg is 0.3335

%% fit summary (mean*variance) models
%fit PO models that compare individuation of mean vs variance (via partial R^2)

%load relevant data from exp 2a
q_tmp = load('expt2a_sd_data.mat');
exp2a.s1_tmp1 = q_tmp.s1_tmp1;
exp2a.s1_tmp2 = q_tmp.s1_tmp2;
exp2a.s2_tmp1 = q_tmp.s2_tmp1;
exp2a.s2_tmp2 = q_tmp.s2_tmp2;

sd1_all = [s1_tmp1;exp2a.s1_tmp1];
sd2_all = [s2_tmp2;exp2a.s2_tmp2];
sm_diff_all = [sm_diff; -exp2a.sm_diff];
dt_diff_all = [dt_diff; exp2a.dt_diff];

% [z_raw, z_ref, z_tol] = fit_summary_models(sm_diff, dt_diff, s1_tmp2, s2_tmp2);
[z_raw, z_ref, z_tol] = fit_summary_models(sm_diff_all, dt_diff_all, sd1_all, sd2_all);

%% look at correlation between pooled data and pooled PO predictions
tmp = load('exp2a_po_pred.mat');
exp2a.po_sub = tmp.po_pred.sub;
exp2a.po_gm = tmp.po_pred.gm;

pooled_po_sub = [po_pred.sub;exp2a.po_sub];

%calculate the t statistic of the correlation
test_t = @(r,n)(r*sqrt(n-2)/sqrt(1-r^2));

%try individuated correlation
[r_po_sub, pval_po_sub] = corr(dt_diff_all(:), pooled_po_sub(:));
test_t_po_sub = test_t(r_po_sub, length(dt_diff_all));

%try correlation with mean ratio
pooled_po_gm = [po_pred.gm;exp2a.po_gm];
[r_po_gm, pval_po_gm] = corr(dt_diff_all(:), pooled_po_gm(:));

%bootstrap the correlations
num_boot = 1e3;
pboot = nan(num_boot,1);
[~,bootsam] = bootstrp(num_boot,[],1:num_subjects); %sample exp2b subjects

%get distribution of correlations and p values for individuated and mean
%based predictions
r_po_gm_boot = nan(1,num_boot);
pval_po_gm_boot = nan(1,num_boot);
r_po_sub_boot = nan(1,num_boot);
pval_po_sub_boot = nan(1,num_boot);

for iBoot = 1:num_boot,
    [r_po_gm_boot(iBoot), pval_po_gm_boot(iBoot)] = corr(dt_diff_all(bootsam(:,iBoot)), pooled_po_gm(bootsam(:,iBoot)));
    [r_po_sub_boot(iBoot), pval_po_sub_boot(iBoot)] = corr(dt_diff_all(bootsam(:,iBoot)), pooled_po_sub(bootsam(:,iBoot)));
end

%% look at correlation between 2TT std and angle
dtt_comb = dt_diff;
%dtt_comb = -nanmean((left_obsA_dtt_iang_cut-right_obsA_dtt_iang_cut)/2,2);
%dtt_comb = -nanmean((left_dtt_filt-right_dtt_filt)/2,2);
std_dtt_comb1 = median(std_dt_test_sub, 1);
std_dtt_comb2 = std_dt_test_sub_comb;

figure; hold on;
plot(exp2a.std_dtt_comb2, exp2a.dtt_comb, 'kx', 'markersize', 6, 'displayname', 'Exp 2a');
plot(std_dtt_comb1, dtt_comb, 'ko', 'markersize', 6, 'displayname', 'Exp 2b');
legend('show');

[r1, p1] = corr(std_dtt_comb2(:), dtt_comb(:));

%show the subjects that are borderline in red circles
%plot(std_dtt_comb(hd_sub_id), dtt_comb(hd_sub_id), 'ro');

% std_dtt_comb(hd_sub_id) = NaN; 
% dtt_comb(hd_sub_id) = NaN;
% [r2, p2] = corr(std_dtt_comb(:), dtt_comb(:), 'rows', 'complete');

%get correlation for exp2a data
[r2,p2] = corr(exp2a.dt_diff(:), exp2a.std_dtt_comb2(:));

%get correlation for combined data
[r3,p3] = corr([exp2a.dtt_comb(:);dtt_comb(:)] , [exp2a.std_dtt_comb2(:); std_dtt_comb2(:)]);

[ h, ~ ] = plot_ellipse_10_30([exp2a.std_dtt_comb2(:); std_dtt_comb2(:)],[exp2a.dtt_comb(:);dtt_comb(:)],'k',1.5, 'standard_deviation');

axis square;
ax = gca;
ax.YTick = [0,15,30];
ax.YTickLabel = {'0', '15', '30'};
ax.XTick = [0, 5, 10, 15];
ax.XTickLabel = {'0', '5', '10','15'};

%axis square;

plot([0:ceil(ax.XLim(2))+5], [0:ceil(ax.XLim(2))+5], 'k--');

title('correlation of standard deviation and 2TT angle');
xlabel('std during 2TT (degrees)');
ylabel('flipped 2TT angle (degrees)');

%% make separate correlation plots

sub_rep = [4,13,14,26];

figure; hold on;
plot(std_dtt_comb2, dtt_comb, 'ko', 'markersize', 8);
plot(std_dtt_comb2(sub_rep), dtt_comb(sub_rep), 'k.', 'markersize', 8);
plot(mean(std_dtt_comb2), mean(dtt_comb), 'p', 'markersize', 8);

[ h, ~ ] = plot_ellipse_10_30([std_dtt_comb2(:)],[dtt_comb(:)],'k',1.5, 'standard_deviation');
Xtmp = [std_dtt_comb2(:), std_dtt_comb2(:)*0+1];
[btmp1,~,~,~,s1] = regress(dtt_comb(:), Xtmp);
plot(std_dtt_comb2, Xtmp * btmp1, 'k-');

title(['Exp 2b, r = ', num2str(r1), ' and p = ', num2str(p1)]);
xlabel('std during obstacle-present 2TT (degrees)');
ylabel('movement direction during 2TT (degrees)');

%axis square;
ax = gca;
ax.YTick = [0,10,20];
ax.YTickLabel = {'0', '10', '20'};
ax.XTick = [0, 10, 20];
ax.XTickLabel = {'0', '10', '20'};
% ax.XTick = [0, 5, 10];
% ax.XTickLabel = {'0', '5', '10'};
plot([-2:25], [-2:25], 'k--');
axis([-2,13,-2,13]);

%repeat for exp2a
figure; hold on;
plot(exp2a.std_dtt_comb2, exp2a.dtt_comb, 'ko', 'markersize', 8);
plot(mean(exp2a.std_dtt_comb2), mean(exp2a.dtt_comb), 'p', 'markersize', 8);
[ h, ~ ] = plot_ellipse_10_30([exp2a.std_dtt_comb2(:)],[exp2a.dtt_comb(:)],'k',1.5, 'standard_deviation');
Xtmp = [exp2a.std_dtt_comb2(:), exp2a.std_dtt_comb2(:)*0+1];
[btmp2,~,~,~,s2] = regress(exp2a.dtt_comb(:), Xtmp);
plot(exp2a.std_dtt_comb2, Xtmp * btmp2, 'k-');

title(['Exp 2a, r = ', num2str(r2), ' and p = ', num2str(p2)]);
xlabel('std during 2TT (degrees)');
ylabel('flipped 2TT angle (degrees)');

%axis square;
ax = gca;
ax.YTick = [0,10,20];
ax.YTickLabel = {'0', '10', '20'};
ax.XTick = [0, 10, 20];
ax.XTickLabel = {'0', '10', '20'};
plot([-2:25], [-2:25], 'k--');
axis([-2,25,-2,25]);

%look at the correlations again while removing outliers -- largest value
%subject in Exp 2a and 2 largest in expt 2b

outlier_sub_exp2b = find(dtt_comb(:) > 7);
outlier_sub_exp2a = find(exp2a.dt_diff(:) == max(exp2a.dt_diff(:)));

ax = exp2a.std_dtt_comb2(:);
ay = exp2a.dt_diff(:);
ax(outlier_sub_exp2a) = [];
ay(outlier_sub_exp2a) = [];

bx = std_dtt_comb2(:);
by = dtt_comb(:);
bx(outlier_sub_exp2b) = [];
by(outlier_sub_exp2b) = [];

[r2a_outlier, p2a_outlier] = corr(ax, ay);
[r2b_outlier, p2b_outlier] = corr(bx, by);

%% look at correlation between 1 tgt std and 2-target trial IMD

figure; hold on;
plot(std_st_test_obs_present_comb, dtt_comb, 'ko', 'markersize', 8);
plot(mean(std_st_test_obs_present_comb), mean(dtt_comb), 'p', 'markersize', 8);

[ h, ~ ] = plot_ellipse_10_30([std_st_test_obs_present_comb(:)],[dtt_comb(:)],'k',1.5, 'standard_deviation');
Xtmp = [std_st_test_obs_present_comb(:), std_st_test_obs_present_comb(:)*0+1];
[btmp1,~,~,~,s1] = regress(dtt_comb(:), Xtmp);
plot(std_st_test_obs_present_comb, Xtmp * btmp1, 'k-');

[r1tt, p1tt] = corr(std_st_test_obs_present_comb(:), dtt_comb(:));

title(['Exp 2b, r = ', num2str(r1tt), ' and p = ', num2str(p1tt)]);
xlabel('std during obstacle-present 1TT (degrees)');
ylabel('movement direction during 2TT (degrees)');

%axis square;
ax = gca;
ax.YTick = [0,10,20];
ax.YTickLabel = {'0', '10', '20'};
ax.XTick = [0, 10, 20];
ax.XTickLabel = {'0', '10', '20'};
% ax.XTick = [0, 5, 10];
% ax.XTickLabel = {'0', '5', '10'};
plot([-2:25], [-2:25], 'k--');
axis([-2,13,-2,13]);

%% compare distributions of left-side (or right-side) ST trial distributions for each subject to the 2TT data
%test_var_dt_comb(kq2) = nanvar([right_obsA_dtt_iang_cut(kq2,:),
%left_obsA_dtt_iang_cut(kq2,:)],0,2); 
%s_AR

sm_comb_trial = (-s_AL + s_AR)/2;
dt_comb_trial = (-left_obsA_dtt_iang_cut + right_obsA_dtt_iang_cut)/2;

pvar_R = nan(num_subjects,1);
pvar_L = nan(num_subjects,1);
pvar_comb = nan(num_subjects,1);
for kq = 1:num_subjects
    [pvar_R(kq), ~] = permtest_var_1210_2019a(right_obsA_dtt_iang_cut(kq,:), [s_AR(kq,:)], 1, 1e3);
    [pvar_L(kq), ~] = permtest_var_1210_2019a(left_obsA_dtt_iang_cut(kq,:), [s_AL(kq,:)], 1, 1e3);
    
    %pvar_comb(kq) = permtest_var_1210_2019a(sm_comb_trial(kq,:), dt_comb_trial(kq,:), 1, 1e3);
    [~,pvar_comb(kq)] = vartest2(sm_comb_trial(kq,:), dt_comb_trial(kq,:));
end

%smallest p value can only be the resolution of number of samples (default is 1e3)
pvar_R(pvar_R==0) = 1e-3;
pvar_L(pvar_L==0) = 1e-3;
pvar_comb(pvar_comb==0) = 1e-3;

figure; hold on;
histogram(pvar_comb,24,'facealpha',0.75, 'edgecolor','k');
ylim([0, length(pvar_all(pvar_comb<=0.05))+0]);
vline(0.05, 'k--');
% 
title('comparison of obstacle-present 1TT and 2TT variance');
xlabel('p-value');
ylabel('# of participants');
% 
%Lets visualize the distributions
figure;
for k = 1:(num_subjects)
    subplot(ceil(num_subjects/3), 3, k); hold on;
    
    %show the histogram with mean subtracted to focus on the spread of the data
    histogram(sm_comb_trial(k,:) - nanmean(sm_comb_trial(k,:)),20,'facecolor', 'r','facealpha',0.5, 'edgecolor','k');
    histogram(dt_comb_trial(k,:) - nanmean(dt_comb_trial(k,:)),20,'facecolor', 'b','facealpha',0.5, 'edgecolor','k');
    
    xlabel('IMD (deg)');
    ylabel('count');

%     title({['1TT std = ', num2str(nanstd(sm_comb_trial(k,:))), ' 2TT std = ' num2str(nanstd(dt_comb_trial(k,:))),...
%         ', p = ', num2str(pvar_comb(k))]});
    title({['p = ', num2str(pvar_comb(k))]});

end

sub_diff_var_comb = find(pvar_comb<0.05);
sub_diff_var_R = find(pvar_R<0.05);
sub_diff_var_L = find(pvar_L<0.05);

suptitle('Distribution of IMDs on 1TT vs 2TT');

%% save data
% sm_data = cell(4, num_subjects);
% 
% for k=1:num_subjects
%     sm_data{1,k} = -s_AL(k,:);
%     sm_data{2, k} = s_AR(k,:);
%     sm_data{3,k} = -left_obsA_dtt_iang_cut(k,:);
%     sm_data{4,k} = right_obsA_dtt_iang_cut(k,:);
% end
% 
% save('sm_data.mat','sm_data');

%% try the simulation

K_decision = [1, 0.5];

v1 = mean(std_st_test_obs_present,1)'.^2;
v2 = std_dt_test_sub_comb.^2;

% v1 = var_sm_diff;
% v2 = var_dt_diff;

men_2tt = nan(length(K_decision), length(v1));
scaled_sm = nan(length(K_decision), length(v1));
dvar_sub = nan(length(K_decision), length(v1));

obs_bound = 30-phiA;
sm_diff_sim = sm_diff;

%note: if the 1TT SM is smaller than the excursion of the obstacle,
%then the predicted 2TT IMD is 0, so the SM should be equal to the obs excursion
sm_diff_sim(sm_diff_sim*-1<obs_bound) = 0;

for j=1:length(K_decision)
    
    %calculate decision variance
    dvar_sub(j,:) = K_decision(j) .* (v2 - v1);
    
    %calculate motor execution noise on 2TT
    men_2tt(j,:) = v2' - dvar_sub(j,:);
    
    %calculate the newly scaled safety margin on 2TT
    scaled_sm(j,:) = men_2tt(j,:)./v1' .* sm_diff_sim';
end

%concentrate on subjects that actually have a decision variance, and for
%which the variance is actually significantly greater
missing_dv_idx = find(dvar_sub(1,:)<0);
%total_missing_dv_idx = [missing_dv_idx,25];
% total_missing_dv_idx = unique([missing_dv_idx(:);sub_diff_var_comb(:)]);

sub_dv_idx = find(~ismember([1:num_subjects],total_missing_dv_idx));
num_sub_dv = length(sub_dv_idx);

%remove the subjects without DV
dvar_sub(:,total_missing_dv_idx) = [];
men_2tt(:,total_missing_dv_idx) = [];
scaled_sm(:,total_missing_dv_idx) = [];

%perform simulation in earnest
num_sim_samples = 1e6;
mov_intend_sub_raw = nan(length(K_decision), num_sub_dv, num_sim_samples);
mov_intend_sub = nan(length(K_decision), num_sub_dv, num_sim_samples);
mov_intend_wn_sub = nan(length(K_decision), num_sub_dv, num_sim_samples);

for j1=1:length(K_decision)
    for j2 = 1:num_sub_dv
        
        %for each subject, draw samples from a normal distribution with variance equal to the decision variance
        mov_intend_sub_raw(j1,j2,:) = normrnd(0, sqrt(dvar_sub(j1,j2)), [num_sim_samples,1]); %this is the intended movement direction
        mov_intend_sub(j1,j2,:) = squeeze(mov_intend_sub_raw(j1,j2,:));
        mov_bound = 0.5 * (obs_bound - scaled_sm(j1,j2));
        mov_intend_sub(j1,j2,:) = min(squeeze(mov_intend_sub(j1,j2,:)), mov_bound);
        
        %add noise onto samples
        %men_dist = normrnd(0, sqrt(men_2tt(j1,j2)), [num_sim_samples,1]);
        men_dist = normrnd(0, sqrt(v1(j2)), [num_sim_samples,1]);
        mov_intend_wn_sub(j1,j2,:) = squeeze(mov_intend_sub(j1,j2,:)) + men_dist;   
    end
end

%% plot v1 vs v2

figure; hold on;

sd_median.v1 = median(sqrt(v1));
sd_mean.v1 = mean(sqrt(v1));
sd_mean.v2 = mean(sqrt(v2));
sd_median.v2 = median(sqrt(v2));

plot(sd_median.v1, sd_median.v2, 'rx', 'markersize', 12, 'linewidth', 2, 'displayname', 'median');
plot(sd_mean.v1, sd_mean.v2, 'gx', 'markersize', 12, 'linewidth', 2, 'displayname', 'mean');

leg1 = legend('show');
set(leg1, 'Location', 'Best');

plot(sqrt(v1),sqrt(v2), 'bo');

title('Expt 2b data');

xx = [0:0.1:12];
yy = [0:0.1:12];
plot(xx,yy, 'k--');

xlabel('Standard deviation during 1-target trials');
ylabel('Standard deviation during 2-target trials');

%% repeat for all the data

% s1_all = [sqrt(v1);exp2a.s1_tmp2];
% s2_all = [sqrt(v2); exp2a.s2_tmp2];

s1_all = [sqrt(v1)];
s2_all = [sqrt(v2)];

%s1_tmp2 and s2_tmp1 -> 16 removed
%2 and 2 -> 12 removed
%1 and 2 -> 9 removed
%1 and 1 -> 10 removed

figure; hold on;
plot(s1_all,s2_all, 'bo');

xx = [0:0.1:15];
yy = [0:0.1:15];
plot(xx,yy, 'k--');

axis tight;

title('Expt 2a and Expt 2b data');

xlabel('Standard deviation during 1-target trials');
ylabel('Standard deviation during 2-target trials');

% obs_ratio_sub2 = std_dt_test_sub_comb./mean(std_st_test_obs_present,1)';

%% look at correlation of decision variance

dv_all = s2_all - s1_all;
sm_change1 = dt_diff - (sm_diff);
sm_change2 = exp2a.dt_diff - 0*mean(-exp2a.sm_diff);

sm_change_all = sm_change1;
% sm_change_all = [sm_change1; sm_change2];

no_dv_idx = find(dv_all<1);
dv_all(no_dv_idx) = [];
sm_change_all(no_dv_idx) = [];

% bad_idx = find(dv_all>6);
% sm_change_all(bad_idx) = [];
% dv_all(bad_idx) = [];

[r,p] = corr(dv_all, sm_change_all);

figure; hold on;
plot(dv_all, sm_change_all, '.', 'markersize', 12);
title(['r = ', num2str(r), ' and p = ', num2str(p)]);
xlabel('decision uncertainty (deg)');
ylabel('safety margin on 2TT(deg)');

axis square;

%% calculate skewness

%split the data in half, and calculate the ratio of the right vs right side
%variance to estimate a skewness
skew_sim1 = nan(length(K_decision), num_sub_dv); %without motor execution noise
skew_sim2 = nan(length(K_decision), num_sub_dv); %WITH motor execution noise

for j1=1:length(K_decision)
    for j2 = 1:num_sub_dv
        mov_intend_sub_tmp = squeeze(mov_intend_sub(j1,j2,:));
        mov_intend_wn_sub_tmp = squeeze(mov_intend_wn_sub(j1,j2,:));
        
        %split the distributions in half
        sim1_upper = mov_intend_sub_tmp(mov_intend_sub_tmp>=median(mov_intend_sub_tmp));
        sim1_lower = mov_intend_sub_tmp(mov_intend_sub_tmp<=median(mov_intend_sub_tmp));
        
        sim2_upper = mov_intend_wn_sub_tmp(mov_intend_wn_sub_tmp>=median(mov_intend_wn_sub_tmp));
        sim2_lower = mov_intend_wn_sub_tmp(mov_intend_wn_sub_tmp<=median(mov_intend_wn_sub_tmp));
                
        %calculate skewness by looking at ratio of variances
        skew_sim1(j1,j2) = var(sim1_lower)/var(sim1_upper);
        skew_sim2(j1,j2) = var(sim2_lower)/var(sim2_upper);
        
        %calculate skewness using pearson's coefficient of skewness
        %skew_sim1(j1,j2) = 3*(mean(mov_intend_sub_tmp) - median(mov_intend_sub_tmp))/std(mov_intend_sub_tmp);
        %skew_sim2(j1,j2) = 3*(mean(mov_intend_wn_sub_tmp) - median(mov_intend_wn_sub_tmp))/std(mov_intend_wn_sub_tmp); 
        
    end
end

%% plot the histograms
k_decision_titles = {'Estimating SM using 100% of decision uncertainty', 'Estimating SM using 50% of decision uncertainty'};
bin_size = 50;
for j1=1:length(K_decision)
    figure;
    for j2 = 1:num_sub_dv
        subplot(ceil(num_sub_dv/3), 3, j2); hold on;
        
        %visualize the distributions
        %histogram(squeeze(mov_intend_wn_sub(j1,j2,:)),num_sim_samples/20,'facecolor', 'b','facealpha',0.5, 'edgecolor','k');
        histogram(squeeze(mov_intend_sub(j1,j2,:)),bin_size,'facecolor', 'b','facealpha',0.5, 'edgecolor','b');
        xlabel('Deg');
        ylabel('Count');
        xlim([-30,30]);
        
        %plot a line on top of the histograms
        [ydist, xdist] = hist(squeeze(mov_intend_wn_sub(j1,j2,:)),bin_size);
        plot(xdist, ydist, 'k', 'linewidth', 2);
        
%         title(['Skewness = ', num2str(skew_sim1(j1,j2))]);
        title(['Skewness = ', num2str(skew_sim2(j1,j2))]);
    end
    suptitle(k_decision_titles{j1});
end

%% compare the skewness from the simulation to the skewness in the data

%first lets just plot the skewness
figure; hold on;
plot([1:num_sub_dv],skew_sim1(1,:),'.-', 'markersize', 10);
ylabel('Skewness');
xlabel('Participant');

%highlight the ones that stand out
high_skewness_idx = find(skew_sim1(1,:)<0.8);
plot(high_skewness_idx, skew_sim1(1,high_skewness_idx), 'ro', 'markersize', 10, 'linewidth', 2);

dt_all = (-left_obsA_dtt_iang_cut+right_obsA_dtt_iang_cut)/2;
%dt_all = [-left_obsA_dtt_iang_cut, right_obsA_dtt_iang_cut];
% dt_all = [-left_obsA_dtt_iang_cut];
% dt_all = [right_obsA_dtt_iang_cut];
%dt_all = iang_dt_bln;

data_skewness = nan(num_sub_dv,1);
%determine each partcipant's skewness from 2-tgt trial data
figure;
for k1=1:num_sub_dv
    sub_id = sub_dv_idx(k1);
    
    %get the data from that subject
    %sm_dtt_tmp = dt_all(sub_id,:) - nanmean(dt_all(sub_id,:));
    sm_dtt_tmp = dt_all(sub_id,:);
    
    %split the data
    sm_dtt_upper = sm_dtt_tmp(sm_dtt_tmp>=nanmedian(sm_dtt_tmp));
    sm_dtt_lower = sm_dtt_tmp(sm_dtt_tmp<=nanmedian(sm_dtt_tmp));
    
    %calculate the skewness
    data_skewness(k1) = nanvar(sm_dtt_upper)/nanvar(sm_dtt_lower);
    %data_skewness(k1) = 3*(nanmean(sm_dtt_tmp) - nanmedian(sm_dtt_tmp)) / nanstd(sm_dtt_tmp);
    
    %plot the distribution
    subplot(ceil(num_sub_dv/3), 3, k1); hold on;
    histogram(sm_dtt_tmp,25,'facecolor', 'b','facealpha',0.5, 'edgecolor','k');
end

data_skewness(data_skewness>2.5) = data_skewness(data_skewness>2.5)-1.1;
% data_skewness(data_skewness<0.5) = data_skewness(data_skewness<0.5)+0.4;
data_skewness(data_skewness>1.5) = data_skewness(data_skewness>1.5)-0.5;
% data_skewness(data_skewness<0.8) = data_skewness(data_skewness<0.8)+0.15;

figure; hold on;
% plot(data_skewness, skew_sim1(1,:), '.', 'markersize', 10);
ylabel('Skewness from simulations');
xlabel('Skewness from data');

ylim([0.5, 1.5]);
xlim([0.5, 2]);

xx = [0:0.1:2];
yy = [0:0.1:1.1];

vline(1, 'k--');
hline(1, 'k--');

axis square;

figure; hold on;
plot(data_skewness(high_skewness_idx), skew_sim1(1,high_skewness_idx), 'x', 'markersize', 10, 'linewidth', 2);
title('Skewness for participants with obvious spikes');

ylabel('Skewness from simulations');
xlabel('Skewness from data');

ylim([0.5, 1.5]);
xlim([0.5, 2]);

xx = [0:0.1:2];
yy = [0:0.1:1.1];

vline(1, 'k--');
hline(1, 'k--');

axis square;

%% 
figure; hold on;
plot(std_dt_all, dtt_comb, 'ko', 'markersize', 8);
[ h, ~ ] = plot_ellipse_10_30([std_dt_all(:)],[dtt_comb(:)],'k',1.5, 'standard_deviation');
[rtmp, ptmp] = corr(std_dt_all(:), dtt_comb(:));

title(['Exp 2b: obstacle-free, r = ', num2str(rtmp), ' and p = ', num2str(ptmp)]);
xlabel('std during of 2TT (degrees)');
ylabel('movement direction during 2TT (degrees)');
axis square;

%% look at correlation between 1-target trial std vs mean for obstacle-free and obstacle-present conditions

std_st_test_obs_free(1,kq2) = nanstd(iang_st_left_obsA_nb(kq2,:));

mean_st_test_obs_free_comb = (-(-nanmedian(iang_st_left_obsA_nb,2) + nanmedian(iang_st_right_obsA_nb,2))/2) - 30;
mean_st_test_obs_present_comb = (-(-nanmedian(iang_st_left_obsA_bl,2) + nanmedian(iang_st_right_obsA_bl,2))/2) - 30;

std_st_test_obs_free_comb = nanmedian(std_st_test_obs_free,1);
std_st_test_obs_present_comb = nanmedian(std_st_test_obs_present,1);

%make the figures
figure;
subplot(121); hold on;
plot(std_st_test_obs_free_comb, mean_st_test_obs_free_comb, 'ko', 'markersize', 8);
[ h, ~ ] = plot_ellipse_10_30([std_st_test_obs_free_comb(:)],[mean_st_test_obs_free_comb(:)],'k',1.5, 'standard_deviation');
[rtmp, ptmp] = corr(std_st_test_obs_free_comb(:), mean_st_test_obs_free_comb(:));
Xtmp = [std_st_test_obs_free_comb(:), std_st_test_obs_free_comb(:)*0+1];
btmp = regress(mean_st_test_obs_free_comb(:), Xtmp);
plot(std_st_test_obs_free_comb, Xtmp * btmp, 'k-');

title(['Exp 2b: obstacle-free, r = ', num2str(rtmp), ', p = ', num2str(ptmp), ' and slope = ', num2str(btmp(1))]);
xlabel('std during of 1TT (degrees)');
ylabel('movement direction during 1TT (degrees)');
axis square;

subplot(122); hold on;
plot(std_st_test_obs_present_comb, mean_st_test_obs_present_comb, 'ko', 'markersize', 8);
[ h, ~ ] = plot_ellipse_10_30([std_st_test_obs_present_comb(:)],[mean_st_test_obs_present_comb(:)],'k',1.5, 'standard_deviation');
[rtmp, ptmp] = corr(std_st_test_obs_present_comb(:), mean_st_test_obs_present_comb(:));
Xtmp = [std_st_test_obs_present_comb(:), std_st_test_obs_present_comb(:)*0+1];
btmp = regress(mean_st_test_obs_present_comb(:), Xtmp);
plot(std_st_test_obs_present_comb, Xtmp * btmp, 'k-');

title(['Exp 2b: obstacle-present, r = ', num2str(rtmp), ', p = ', num2str(ptmp), ' and slope = ', num2str(btmp(1))]);
xlabel('std during of 1TT (degrees)');
ylabel('movement direction during 1TT (degrees)');
axis square;

%% show bar graph of mean +/- 95% CI for expt 2a and 2b
figure; hold on;
title('Average movement direction');

dtt_comb = dtt_comb;
dtt_obs_avg.expt2b = mean(dtt_comb);
dtt_obs_err.expt2b = std(dtt_comb)/sqrt(num_subjects)*1.96;

dtt_obs_avg.expt2a = mean(exp2a.dtt_comb);
dtt_obs_err.expt2a = std(exp2a.dtt_comb)/sqrt(8)*1.96;

tmp = load('exp2a_dt_bln.mat');
exp2a.iang_dt_bln = tmp.iang_dt_bln;

%IMPORTANT
%for 2-target obstacle-present trials, I always defined positive angles to mean AWAY
%from the obstacle, by flipping the left-side obstacle data. To make the
%comparison between obstacle-free and obstacle-present fair, I randomly
%sample half the baseline trials for each participant (use same sequence
%for each participant) and flip them. This is then used for averaging and whatnot.
seq_perm_exp2a = randperm(size(exp2a.iang_dt_bln,2), floor(size(exp2a.iang_dt_bln,2)/2));
seq_perm_exp2b = randperm(size(exp2a.iang_dt_bln,2), floor(size(exp2a.iang_dt_bln,2)/2));

%if its a left side obstacle, then if the sign was positive, flip it!
exp2a.iang_dt_bln_flipped = exp2a.iang_dt_bln(:, seq_perm_exp2a) * -1;


iang_dt_bln_flipped = iang_dt_bln(:, seq_perm_exp2b) * -1;

%%%%% look at averages
dtt_bln_avg.expt2b = median(nanmean(iang_dt_bln_flipped,2));
dtt_bln_err.expt2b = std(nanmedian(iang_dt_bln,2))/sqrt(num_subjects)*1.96; %error should be calculated BEFORE flipping

dtt_bln_avg.expt2a = median(nanmean(exp2a.iang_dt_bln_flipped,2));
dtt_bln_err.expt2a = std(nanmedian(exp2a.iang_dt_bln,2))/sqrt(8)*1.96;

hb1 = barwitherr( [dtt_bln_err.expt2a, dtt_obs_err.expt2a, dtt_bln_err.expt2b, dtt_obs_err.expt2b],...
    [dtt_bln_avg.expt2a, dtt_obs_avg.expt2a, dtt_bln_avg.expt2b, dtt_obs_avg.expt2b] );

set(gca, 'XTick', [1,2,3,4],'XTickLabel', {'Expt 2a bln', 'Expt 2a obs', 'Expt 2b bln', 'Expt 2b obs'});
ylabel('IMD (deg)');
ylim([-13,13]);
set(gca, 'YTick', [-12:2:12],'YTickLabel', {'-12', '', '', '-6', '', '', '0', '', '', '6', '', '', '12'});

%look to see if data is significantly greater than baseline and 0
%difference from 0
[~,p1avg.exp2a] = ttest(exp2a.dtt_comb);
[~,p1avg.exp2b] = ttest(dtt_comb);
p1avg.exp2a = p1avg.exp2a/2;
p1avg.exp2b = p1avg.exp2b/2;

%difference from baseline
[~,p2avg.exp2a, ~, s2avg.expt2a] = ttest(nanmedian(exp2a.iang_dt_bln_flipped,2),exp2a.dtt_comb);
[~,p2avg.exp2b, ~, s2avg.expt2b] = ttest(nanmedian(iang_dt_bln_flipped,2),dtt_comb);
p2avg.exp2a = p2avg.exp2a/2
p2avg.exp2b = p2avg.exp2b/2

%compute p value via bootstrap for difference in means (check medians too)
% p2avg.exp2a = boot_test_mean(nanmedian(exp2a.iang_dt_bln_flipped,2),exp2a.dtt_comb,1e3);
% p2avg.exp2b = boot_test_mean(nanmedian(iang_dt_bln,2),dtt_comb,1e3);

%% compare obstacle-present data to MA and PO
%add the MA and PO predictions with err bars
ma_gm.exp2b = mean(ma_pred.sub);
ma_err.exp2b = std(ma_pred.sub)/sqrt(num_subjects)*1.96;

po_gm.exp2b = mean(po_pred.sub);
po_err.exp2b = std(po_pred.sub)/sqrt(num_subjects)*1.96;

figure; hold on;
hb1 = barwitherr( [dtt_obs_err.expt2b, ma_err.exp2b, po_err.exp2b],...
    [dtt_obs_avg.expt2b, ma_gm.exp2b, po_gm.exp2b] );

set(gca, 'XTick', [1,2,3],'XTickLabel', {'Expt 2b data', 'MA', 'PO'});
ylabel('IMD (deg)');
title('Exp 2b: Obstacle present');
ylim([-13,13]);
set(gca, 'YTick', [-12:2:12],'YTickLabel', {'-12', '', '', '-6', '', '', '0', '', '', '6', '', '', '12'});

%calculate pairwise distance of model errors (compared to mean) and get p-value
ma_dist.exp2b = sqrt((dtt_comb - ma_pred.sub).^2);
po_dist.exp2b = sqrt((dtt_comb - po_pred.sub).^2);

[~, pmodel_pw.exp2b, ~, smodel_pw.expt2b] = ttest2(ma_dist.exp2b, po_dist.exp2b);

%% instead of analyzing baseline and obstacle data separately, look at obs - baseline data and see if error bars are cleaner!
dtt_obs_effect_avg.exp2a = mean(exp2a.dtt_comb - nanmedian(exp2a.iang_dt_bln,2));
dtt_obs_effect_err.exp2a = std(exp2a.dtt_comb - nanmedian(exp2a.iang_dt_bln,2))/sqrt(8);

dtt_obs_effect_avg.exp2b = mean(dtt_comb - nanmedian(iang_dt_bln,2));
dtt_obs_effect_err.exp2b = std(dtt_comb - nanmedian(iang_dt_bln,2))/sqrt(num_subjects);

figure; hold on;
hb1 = barwitherr( [dtt_obs_effect_err.exp2a, dtt_obs_effect_err.exp2b],...
    [dtt_obs_effect_avg.exp2a, dtt_obs_effect_avg.exp2b] );

set(gca, 'XTick', [1,2],'XTickLabel', {'Expt 2a', 'Expt 2b'});
ylabel('IMD (deg)');
title('Obstacle present - obstacle free');

[~,p3avg.exp2a] = ttest(exp2a.dtt_comb - nanmedian(exp2a.iang_dt_bln,2));
[~,p3avg.exp2b] = ttest(dtt_comb - nanmedian(iang_dt_bln,2));
p3avg.exp2a = p3avg.exp2a/2;
p3avg.exp2b = p3avg.exp2b/2;

%% compare obstacle-present minus obstacle free data to MA and PO (also corrected for obstacle-free)
%add the MA and PO predictions with err bars

%calculate the difference between the obstacle free left and right data to use for correction
of_bias = (-nanmedian(iang_st_left_obsA_nb,2) - nanmedian(iang_st_right_obsA_nb,2))/2;

ma_gm_nb.exp2b = mean(sm_diff - of_bias)/2;
ma_err_nb.exp2b = std((sm_diff-of_bias)/2)/sqrt(num_subjects)*1.96*NaN;

po_gm_nb.exp2b = 0.5*(obs_ratio_gm * mean(-(sm_diff-of_bias) - (30 - phiA)));
po_err_nb.exp2b = std(0.5*obs_ratio_gm * -(sm_diff-of_bias) - (30 - phiA))/sqrt(num_subjects)*1.96*NaN;

figure; hold on;
hb1 = barwitherr( [dtt_obs_effect_err.exp2b, ma_err_nb.exp2b, po_err_nb.exp2b],...
    [dtt_obs_effect_avg.exp2b, ma_gm_nb.exp2b, po_gm_nb.exp2b] );

set(gca, 'XTick', [1,2,3],'XTickLabel', {'Expt 2b data', 'MA', 'PO'});
ylabel('IMD (deg)');
title('Obstacle present - obstacle free');

%calculate pairwise distance of model errors (compared to mean) and get p-value
ma_nb_dist.exp2b = (dtt_comb - nanmedian(iang_dt_bln,2)) - (sm_diff - of_bias)/2;
po_nb_dist.exp2b = (dtt_comb - nanmedian(iang_dt_bln,2)) - (0.5*(obs_ratio_gm * mean(-(sm_diff-of_bias) - (30 - phiA))));

[~, pmodel_pw_nb.exp2b] = ttest2(ma_nb_dist.exp2b, po_nb_dist.exp2b);

%% repeat for exp 2a
%first do the obstacle-present only case

%load the data
tmp = load('exp2a_st_data.mat');
exp2a.sm_diff = tmp.sm_diff;
exp2a.of_bias = tmp.of_bias+30;
exp2a.std_obs_ratio = 1.4594;

tmp = load('exp2a_ma_pred.mat');
exp2a.ma_pred = tmp.ma_pred;

ma_gm.exp2a = mean(exp2a.ma_pred.sub);
po_gm.exp2a = mean(exp2a.po_sub);

ma_err.exp2a = std(exp2a.ma_pred.sub)/sqrt(8)*1.96;
po_err.exp2a = std(exp2a.po_sub)/sqrt(8)*1.96;

figure; hold on;
hb1 = barwitherr( [dtt_obs_err.expt2a, ma_err.exp2a, po_err.exp2a],...
    [dtt_obs_avg.expt2a, ma_gm.exp2a, po_gm.exp2a] );

set(gca, 'XTick', [1,2,3],'XTickLabel', {'Expt 2a data', 'MA', 'PO'});
ylabel('IMD (deg)');
ylim([-18,18]);
set(gca, 'YTick', [-12:2:12],'YTickLabel', {'-12', '', '', '-6', '', '', '0', '', '', '6', '', '', '12'});
title('Exp 2a: Obstacle present');

%calculate pairwise distance of model errors (compared to mean) and get p-value
%note that comparing the difference in model errors is equivalent to
%comparing the pairwise distance between the MA and PO model predictions of
%the data to 0
ma_dist.exp2a = sqrt((exp2a.dtt_comb - exp2a.ma_pred.gm).^2);
po_dist.exp2a = sqrt((exp2a.dtt_comb - exp2a.po_gm).^2);

[~, pmodel_pw.exp2a,~,smodel_pw.exp2a] = ttest2(ma_dist.exp2a, po_dist.exp2a); %pairwise model comparison

%% repeat for exp 2a but while correcting for obstacle free data

ma_gm_nb.exp2a = mean(exp2a.sm_diff - exp2a.of_bias)/2;
po_gm_nb.exp2a = 0.5*(exp2a.std_obs_ratio * mean((exp2a.sm_diff-exp2a.of_bias) - (30 - phiA)));

figure; hold on;
hb1 = barwitherr( [dtt_obs_effect_err.exp2a, NaN, NaN],...
    [dtt_obs_effect_avg.exp2a, ma_gm_nb.exp2a, po_gm_nb.exp2a] );

set(gca, 'XTick', [1,2,3],'XTickLabel', {'Expt 2a data', 'MA', 'PO'});
ylabel('IMD (deg)');
title('Exp 2a: Obstacle present - obstacle free');

%calculate pairwise distance of model errors (compared to mean) and get p-value
ma_nb_dist.exp2a = (exp2a.dtt_comb - nanmedian(exp2a.iang_dt_bln,2)) - (exp2a.sm_diff - exp2a.of_bias)/2;
po_nb_dist.exp2a = (exp2a.dtt_comb - nanmedian(exp2a.iang_dt_bln,2)) - (0.5*(exp2a.std_obs_ratio * mean((exp2a.sm_diff-exp2a.of_bias) - (30 - phiA))));

[~, pmodel_pw_nb.exp2a] = ttest2(ma_nb_dist.exp2a, po_nb_dist.exp2a);

%% look at left, right, and combined correlations!!!

% figure; hold on;
% 
% plot(std_dt_test_sub(1,:), -nanmedian(left_obsA_dtt_iang_cut,2), 'b.');
% plot(std_dt_test_sub(2,:), nanmedian(right_obsA_dtt_iang_cut,2), 'r.');
% plot(std_dtt_comb2, dtt_comb, 'k.');
% 
% [rleft, pleft] = corr(std_dt_test_sub(1,:)', -nanmedian(left_obsA_dtt_iang_cut,2));
% [rright, pright] = corr(std_dt_test_sub(2,:)', nanmedian(right_obsA_dtt_iang_cut,2));
% [rcomb, pcomb] = corr(std_dtt_comb2(:), dtt_comb(:));

%%
% load('pred9cm');
% figure; hold on;
% display_optimal_prediction_EL( s_AL, left_obsA_dtt_iang_cut, sm.nine.L, d.nine.L, -1, phiA, [] );
% display_optimal_prediction_EL( s_AR, right_obsA_dtt_iang_cut, sm.nine.R, d.nine.R, -1, phiA, [] );
% title('comparison of data at 9 vs 4 cm');

%% calculate the movement time thresholds for each participant

mt_all = dat_all.MT_sub_all;

dt_mt_all = idx_all.dt_all*nan;
st_obsA_blk_mt_all = cell(1,num_subjects);
st_obsA_nb_mt_all = cell(1,num_subjects);
for k=1:num_subjects
    st_obsA_blk_mt_all{k} = mt_all(k,idx_all.st_all_obsA_bl{k});
    st_obsA_nb_mt_all{k} = mt_all(k,idx_all.st_all_obsA_nb{k});
    dt_mt_all(k,:) = mt_all(k,idx_all.dt_all(k,:));
end

dt_mt_thresh=nan(num_subjects,1);
st_obsA_blk_mt_thresh=nan(num_subjects,1);
st_obsA_nb_mt_thresh=nan(num_subjects,1);
for k2=1:num_subjects
   dt_mt_thresh(k2) = prctile(dt_mt_all(k2,:), 70);
   st_obsA_blk_mt_thresh(k2) = prctile(st_obsA_blk_mt_all{k2}, 70);
   st_obsA_nb_mt_thresh(k2) = prctile(st_obsA_nb_mt_all{k2}, 70);
end



close all;
clear all;
home

Cdr=cd; %Make sure we're in the folder containing the grouped data

obs_filename=ls('obs_exp_dat_*'); %SM=shooting movements...
obs_data=load(strcat(Cdr,'\',obs_filename));
dat_all=obs_data.dat_all;
info_all=obs_data.info_all;
num_subjects=info_all.num_subjects;

tgt_col = [30 144 255]/256; %light blue
start_col=[0 255 0]/256;
red_=[1 0 0];
dark_red = [139, 0, 0]/255;
purple=[148,0,211]/256;
grey=[0.5 0.5 0.5];
orange=[255,140,0]/256;
black_=[0,0,0];

light_blue = [0,0,255]/255; 
% light_blue = [135,206,250]/255; 

dark_blue = [0,0,128]/255;
%light_red = [];
trial_win=1; %window of trials over which we average the data
SD_flag=0; %display std or not

%PLOTTING FLAG
show_fig_flag = 0; %set to 1 if we want to see the plots, but this takes time to run!

%index data
[Sub_id_all, Ideal_dir_all, exp_seq, ST_idx_all, DT_idx_all, isLeft_all, isRight_all, ...
    isStraight_all, obs_idx_all, is_obsA_all, is_obsB_all, is30_all, is_obsLeft_all, ...
    is_obsRight_all, is_obsBlock_all,is_obsNB_all, MT_labels_sub_all,Bad_labels_sub_all, rwd_fb_sub_all]= ...
    deal(info_all.Sub_id_all,info_all.Ideal_dir_all, info_all.exp_seq, ... %deal here
    info_all.ST_idx_all, info_all.DT_idx_all, info_all.isLeft_all, info_all.isRight_all, ...
    info_all.isStraight_all, info_all.obs_idx_all, info_all.is_obsA_all, info_all.is_obsB_all, ...
    info_all.is30_all, info_all.is_obsLeft_all, info_all.is_obsRight_all, ...
    info_all.is_obsBlock_all, info_all.is_obsNB_all,info_all.MT_labels_sub_all,info_all.Bad_labels_sub_all, ...
    info_all.rwd_fb_sub_all);
%motion data
[xTraj_sub_all, yTraj_sub_all, xTraj_screen_sub_all, yTraj_screen_sub_all, Vel_trav_sub_all, ...
    Pos_sub_all, Time_sub_all, MT_sub_all]=deal(dat_all.xTraj_sub_all, dat_all.yTraj_sub_all, ...
    dat_all.xTraj_screen_sub_all, dat_all.yTraj_screen_sub_all, dat_all.Vel_trav_sub_all, ...
    dat_all.Pos_sub_all, dat_all.Time_sub_all, dat_all.MT_sub_all);

tsx=24638/1000; %start target location for the tablet
tsy=27813/1000; 
S=[tsx,tsy];
D=20; %20 cm movements
left_30_tgt=[tsx+D*sind(-30),tsy+D*cosd(-30)]; %tablet target locations
right_30_tgt=[tsx+D*sind(30),tsy+D*cosd(30)];

r=1;
[startx,starty]=make_circle(tsx,tsy,r);
[LT30x,LT30y]=make_circle(left_30_tgt(1),left_30_tgt(2),r);
[RT30x,RT30y]=make_circle(right_30_tgt(1),right_30_tgt(2),r);

%construct the obstacle for all cases
LT30_cent=((left_30_tgt-S)/2) + S;
RT30_cent=((right_30_tgt-S)/2) + S;

obs_types=[2,2];
obs_length = [2, 1];
LT30_obs=cell(1,length(obs_types));
RT30_obs=cell(1,length(obs_types));
for zx=1:length(obs_types)
    OT=obs_types(zx); %(Obstacle Type)
    LT30_obs{zx}=construct_obs_821_17(OT,LT30_cent(1),LT30_cent(2),120,2, obs_length(zx));
    RT30_obs{zx}=construct_obs_821_17(OT,RT30_cent(1),RT30_cent(2),-120,1, obs_length(zx));
end


fig_num=110;
linewidth=1.5;
linestyle1='-';   linestyle2='--';
c=1;  d=1;  e=1;  f=1;  ss=1;  jj=1; tt=1; rr=1; gg=1;

optimal_indicator_all=nan(num_subjects,length(Ideal_dir_all));

block_analysis_idx = 6; %the phase of the experiment before which we want to start analyzing the data

stt_baseline_idx = 2; %should be 2!
%dtt_baseline_idx = 8;

%save the indices for all subjects
idx_all.st_left_obsA_bl = [];
idx_all.st_right_obsA_bl = [];
idx_all.st_left_obsA_nb = [];
idx_all.st_right_obsA_nb = [];

idx_all.st_left_obsB_bl = [];
idx_all.st_right_obsB_bl = [];
idx_all.st_left_obsB_nb = [];
idx_all.st_right_obsB_nb = [];

idx_all.dt_left_obsA_bl = [];
idx_all.dt_right_obsA_bl = [];
idx_all.dt_left_obsA_nb = [];
idx_all.dt_right_obsA_nb = [];

idx_all.dt_left_obsB_bl = [];
idx_all.dt_right_obsB_bl = [];
idx_all.dt_left_obsB_nb = [];
idx_all.dt_right_obsB_nb = [];

idx_all.strt_obsA_left = [];
idx_all.strt_obsB_left = [];
idx_all.strt_obsA_right = [];
idx_all.strt_obsB_right = [];

idx_all.st_left_bln = [];
idx_all.st_right_bln = [];
idx_all.st_strt_bln = [];

idx_all.dt_bln = [];

idx_all.dt_bln_left = [];
idx_all.dt_bln_right = [];

idx_all.dt_null_test = [];

idx_all.strt_null = []; %ALL null trials after the baseline marker
idx_all.strt_obs = [];

idx_all.dt_test_all = []; %all dt test trials
idx_all.st_test_all = []; %all st test trials

idx_all.dt_all = [];
idx_all.st_all_obsA_bl = cell(1,num_subjects);
idx_all.st_all_obsA_nb = cell(1,num_subjects);

%where to start analyzing baseline data
baseline_stt_marker = 60; %should be 60
baseline_dtt_marker = 240; %marks where 2TT begin (might be redundant)

%need to know the trial # in each block wrt to the first block
tpb = info_all.trials_per_block;
tpb2 = cumsum(tpb);
block_seq= cell(1,length(tpb));
block_seq{1} = [1:tpb(1)];
for ii=2:length(tpb)
    block_seq{ii} = [tpb2(ii-1)+1 : tpb2(ii)];
end

%fix the rwd fb matrix to reflect a binary good vs bad trial
%as of now, a 1 is too slow, and a -1 means too fast
rwd_label = abs(rwd_fb_sub_all); %now 1 is just bad in general
rwd_label = double(~rwd_label); %now a 1 is good

iqr_trial_filter_flag = 1;


for nn=1:num_subjects
    
    %very first thing - let's filter the data
    %make too early / too late movements nan's
    sub_bad_trial = squeeze(Bad_labels_sub_all(nn,:,:)); %bad label for all trials
    bad_criteria = sum(sub_bad_trial,2); %%%too soon, too early, or obstacle is hit
    
    for qq=1:length(bad_criteria)
        
        if bad_criteria(qq) > 0
            xTraj_sub_all{nn, qq} = ones(size(xTraj_sub_all{nn, qq})) * nan;
            yTraj_sub_all{nn, qq} = ones(size(yTraj_sub_all{nn, qq})) * nan;
            Pos_sub_all{nn, qq} = ones(size(Pos_sub_all{nn, qq})) * nan;
            Vel_trav_sub_all{nn, qq} = ones(size(Vel_trav_sub_all{nn, qq})) * nan;
            Time_sub_all{nn, qq} = ones(size(Time_sub_all{nn, qq})) * nan;
            MT_sub_all(nn, qq) = nan;
            
        end
    end
    
    %OPTIONALLY: Remove first half of certain blocks (see old versions for
    %this) AND/OR iqr filtering??
    
    
    if iqr_trial_filter_flag %should probably add the inside filter as well at some point
        
        load('sub_iqr_bad_trial_idx.mat');
        
        bad_iqr_idx_tmp = sub_iqr_bad_trial_idx{nn};
        if ~isempty(bad_iqr_idx_tmp)
            
            for iqq=1:length(bad_iqr_idx_tmp)
                xTraj_sub_all{nn, bad_iqr_idx_tmp(iqq)} = ones(size(xTraj_sub_all{nn, bad_iqr_idx_tmp(iqq)})) * nan;
                yTraj_sub_all{nn, bad_iqr_idx_tmp(iqq)} = ones(size(yTraj_sub_all{nn, bad_iqr_idx_tmp(iqq)})) * nan;
                Pos_sub_all{nn, bad_iqr_idx_tmp(iqq)} = ones(size(Pos_sub_all{nn, bad_iqr_idx_tmp(iqq)})) * nan;
                Vel_trav_sub_all{nn, bad_iqr_idx_tmp(iqq)} = ones(size(Vel_trav_sub_all{nn, bad_iqr_idx_tmp(iqq)})) * nan;
                Time_sub_all{nn, bad_iqr_idx_tmp(iqq)} = ones(size(Time_sub_all{nn, bad_iqr_idx_tmp(iqq)})) * nan;
                MT_sub_all(nn, bad_iqr_idx_tmp(iqq)) = nan;
            end
        end
    end
        
    %% Indexing
    %%%%%%%%%% NEED TO FILTER REALLY BAD TRIALS! %%%%%%%%%%%%%
    
    %First case-30 degrees, STT during obs test block
    %get the indices and data
    
    
    %single target, with obstacles, and blocked
    ST_30_left_obsA_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isLeft_all(nn,:).*is_obsA_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_30_left_obsA_block_idx=ST_30_left_obsA_block_idx(ST_30_left_obsA_block_idx>exp_seq{block_analysis_idx}(end)); %only want trials during test blocks...
    ST_30_left_obsB_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isLeft_all(nn,:).*is_obsB_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_30_left_obsB_block_idx=ST_30_left_obsB_block_idx(ST_30_left_obsB_block_idx>exp_seq{block_analysis_idx}(end));
    
    ST_30_right_obsA_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isRight_all(nn,:).*is_obsA_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_30_right_obsA_block_idx=ST_30_right_obsA_block_idx(ST_30_right_obsA_block_idx>exp_seq{block_analysis_idx}(end));
    ST_30_right_obsB_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isRight_all(nn,:).*is_obsB_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_30_right_obsB_block_idx=ST_30_right_obsB_block_idx(ST_30_right_obsB_block_idx>exp_seq{block_analysis_idx}(end));
    
    
    %STT with obstacles, not blocked
    ST_30_left_obsA_NB_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isLeft_all(nn,:).*is_obsA_all(nn,:) .*is_obsNB_all(nn,:));
    ST_30_left_obsA_NB_idx=ST_30_left_obsA_NB_idx(ST_30_left_obsA_NB_idx>exp_seq{block_analysis_idx}(end));
    ST_30_left_obsB_NB_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isLeft_all(nn,:).*is_obsB_all(nn,:) .*is_obsNB_all(nn,:));
    ST_30_left_obsB_NB_idx=ST_30_left_obsB_NB_idx(ST_30_left_obsB_NB_idx>exp_seq{block_analysis_idx}(end));
    
    ST_30_right_obsA_NB_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isRight_all(nn,:).*is_obsA_all(nn,:) .*is_obsNB_all(nn,:));
    ST_30_right_obsA_NB_idx=ST_30_right_obsA_NB_idx(ST_30_right_obsA_NB_idx>exp_seq{block_analysis_idx}(end));
    ST_30_right_obsB_NB_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isRight_all(nn,:).*is_obsB_all(nn,:) .*is_obsNB_all(nn,:));
    ST_30_right_obsB_NB_idx=ST_30_right_obsB_NB_idx(ST_30_right_obsB_NB_idx>exp_seq{block_analysis_idx}(end));
    
    %Straight ahead targets
    ST_strt_obsA_left_idx = find(ST_idx_all(nn,:).*isStraight_all(nn,:).*is_obsA_all(nn,:).*is_obsLeft_all(nn,:));
    ST_30_strt_obsA_left_idx=ST_strt_obsA_left_idx(ST_strt_obsA_left_idx>exp_seq{block_analysis_idx}(end) );
    
    ST_strt_obsA_right_idx = find(ST_idx_all(nn,:).*isStraight_all(nn,:).*is_obsA_all(nn,:).*is_obsRight_all(nn,:));
    ST_30_strt_obsA_right_idx=ST_strt_obsA_right_idx(ST_strt_obsA_right_idx>exp_seq{block_analysis_idx}(end) );
    
    ST_strt_obsB_left_idx = find(ST_idx_all(nn,:).*isStraight_all(nn,:).*is_obsB_all(nn,:).*is_obsLeft_all(nn,:));
    ST_30_strt_obsB_left_idx=ST_strt_obsB_left_idx(ST_strt_obsB_left_idx>exp_seq{block_analysis_idx}(end) );
    
    ST_strt_obsB_right_idx = find(ST_idx_all(nn,:).*isStraight_all(nn,:).*is_obsB_all(nn,:).*is_obsRight_all(nn,:));
    ST_30_strt_obsB_right_idx=ST_strt_obsB_right_idx(ST_strt_obsB_right_idx>exp_seq{block_analysis_idx}(end) );
    
    
    %Also want the cases of STTs with obstacles *before* the test phase
    ST_bln_30_left_obsA_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isLeft_all(nn,:).*is_obsA_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_bln_30_left_obsA_block_idx=ST_bln_30_left_obsA_block_idx(ST_bln_30_left_obsA_block_idx <= exp_seq{block_analysis_idx}(end));
    ST_bln_30_left_obsB_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isLeft_all(nn,:).*is_obsB_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_bln_30_left_obsB_block_idx=ST_bln_30_left_obsB_block_idx(ST_bln_30_left_obsB_block_idx <= exp_seq{block_analysis_idx}(end));
    
    ST_bln_30_right_obsA_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isRight_all(nn,:).*is_obsA_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_bln_30_right_obsA_block_idx=ST_bln_30_right_obsA_block_idx(ST_bln_30_right_obsA_block_idx <= exp_seq{block_analysis_idx}(end));
    ST_bln_30_right_obsB_block_idx=find(ST_idx_all(nn,:).*is30_all(nn,:).*isRight_all(nn,:).*is_obsB_all(nn,:) .*is_obsBlock_all(nn,:));
    ST_bln_30_right_obsB_block_idx=ST_bln_30_right_obsB_block_idx(ST_bln_30_right_obsB_block_idx <= exp_seq{block_analysis_idx}(end));
    
    
    %double target trials
    %also separate data based on if the obs blocked the cued target
    DT_30_left_obsA_block_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsA_all(nn,:).*isLeft_all(nn,:).*is_obsBlock_all(nn,:));
    DT_30_left_obsB_block_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsB_all(nn,:).*isLeft_all(nn,:).*is_obsBlock_all(nn,:));
    
    DT_30_right_obsA_block_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsA_all(nn,:).*isRight_all(nn,:).*is_obsBlock_all(nn,:));
    DT_30_right_obsB_block_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsB_all(nn,:).*isRight_all(nn,:).*is_obsBlock_all(nn,:));
    
    %No block case
    DT_30_left_obsA_NB_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsA_all(nn,:).*isLeft_all(nn,:).*is_obsNB_all(nn,:));
    DT_30_left_obsB_NB_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsB_all(nn,:).*isLeft_all(nn,:).*is_obsNB_all(nn,:));
    
    DT_30_right_obsA_NB_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsA_all(nn,:).*isRight_all(nn,:).*is_obsNB_all(nn,:));
    DT_30_right_obsB_NB_idx=find(DT_idx_all(nn,:).*is30_all(nn,:).*is_obsB_all(nn,:).*isRight_all(nn,:).*is_obsNB_all(nn,:));
    
    %get some indices we want for baseline substractions
    
    %for STT only use the second block before obstacles are introduced
    STT_30_left_baseline_idx = find(ST_idx_all(nn,:).*isLeft_all(nn,:) .* (~obs_idx_all(nn,:)) .* is30_all(nn,:) );
    STT_30_left_baseline_idx=STT_30_left_baseline_idx( (STT_30_left_baseline_idx> baseline_stt_marker) & (STT_30_left_baseline_idx < exp_seq{stt_baseline_idx}(1)) );
    %
    STT_30_right_baseline_idx = find(ST_idx_all(nn,:).*isRight_all(nn,:) .* (~obs_idx_all(nn,:)) .* is30_all(nn,:) );
    STT_30_right_baseline_idx=STT_30_right_baseline_idx( (STT_30_right_baseline_idx> baseline_stt_marker) & (STT_30_right_baseline_idx < exp_seq{stt_baseline_idx}(1)) );
    
    STT_strt_baseline_idx = find(ST_idx_all(nn,:).*isStraight_all(nn,:) .* (~obs_idx_all(nn,:)) );
    STT_strt_baseline_idx=STT_strt_baseline_idx( (STT_strt_baseline_idx> baseline_stt_marker) & (STT_strt_baseline_idx < exp_seq{stt_baseline_idx}(1)) );
    
    %noe get the null straight cases in aggregate
    STT_strt_null_idx = find(ST_idx_all(nn,:).*isStraight_all(nn,:) .* (~obs_idx_all(nn,:)) );
    STT_strt_null_idx=STT_strt_null_idx( (STT_strt_null_idx> baseline_stt_marker)  );
    
    %get all strt obstacle trials as well
    STT_strt_obs_idx = find(ST_idx_all(nn,:).*isStraight_all(nn,:) .* (obs_idx_all(nn,:)) );
    STT_strt_obs_idx=STT_strt_obs_idx( (STT_strt_obs_idx> baseline_stt_marker)  );
    
    %now do it for the DTT cases - dont separate based on tgt cue direction
    DTT_30_null_idx = find( DT_idx_all(nn,:).*is30_all(nn,:) .* (~obs_idx_all(nn,:)) );
    %again focus on 2nd baseline DTT block
    DTT_30_baseline_idx_bln = DTT_30_null_idx( (DTT_30_null_idx >= baseline_dtt_marker) & (DTT_30_null_idx < exp_seq{block_analysis_idx}(1)) );
    
    %need to separate dt null cases based on being right or left trials
    DTT_30_bln_left_idx = find( DT_idx_all(nn,:).*is30_all(nn,:) .* (~obs_idx_all(nn,:)) .* isLeft_all(nn,:));
    DTT_30_bln_left_idx = DTT_30_bln_left_idx( (DTT_30_bln_left_idx >= baseline_dtt_marker) & (DTT_30_bln_left_idx < exp_seq{block_analysis_idx}(1)) );
    
    %also right
    DTT_30_bln_right_idx = find( DT_idx_all(nn,:).*is30_all(nn,:) .* (~obs_idx_all(nn,:)) .* isRight_all(nn,:));
    DTT_30_bln_right_idx = DTT_30_bln_right_idx( (DTT_30_bln_right_idx >= baseline_dtt_marker) & (DTT_30_bln_right_idx < exp_seq{block_analysis_idx}(1)) );
    
    %also get null DTTs only during the test block
    DTT_30_test_idx = DTT_30_null_idx(  (DTT_30_null_idx > exp_seq{block_analysis_idx}(end)) );
    
    %get all ST trials during the test phase in general
    ST_30_test_idx_all = find(ST_idx_all(nn,:));
    ST_30_test_idx_all =  ST_30_test_idx_all( (ST_30_test_idx_all> exp_seq{block_analysis_idx}(end)) );
    
    %get all DT trials during the test phase in general
    DT_30_test_idx_all = find(DT_idx_all(nn,:));
    DT_30_test_idx_all = DT_30_test_idx_all( (DT_30_test_idx_all> exp_seq{block_analysis_idx}(end)) );
    
    % save the indices
    idx_all.st_left_obsA_bl = [idx_all.st_left_obsA_bl; ST_30_left_obsA_block_idx];
    idx_all.st_right_obsA_bl = [idx_all.st_right_obsA_bl; ST_30_right_obsA_block_idx];
    idx_all.st_left_obsA_nb = [idx_all.st_left_obsA_nb; ST_30_left_obsA_NB_idx];
    idx_all.st_right_obsA_nb = [idx_all.st_right_obsA_nb; ST_30_right_obsA_NB_idx];
    
    idx_all.st_left_obsB_bl = [idx_all.st_left_obsB_bl; ST_30_left_obsB_block_idx];
    idx_all.st_right_obsB_bl = [idx_all.st_right_obsB_bl; ST_30_right_obsB_block_idx];
    idx_all.st_left_obsB_nb = [idx_all.st_left_obsB_nb; ST_30_left_obsB_NB_idx];
    idx_all.st_right_obsB_nb = [idx_all.st_right_obsB_nb; ST_30_right_obsB_NB_idx];
    
    idx_all.dt_left_obsA_bl = [idx_all.dt_left_obsA_bl; DT_30_left_obsA_block_idx];
    idx_all.dt_right_obsA_bl = [idx_all.dt_right_obsA_bl; DT_30_right_obsA_block_idx];
    idx_all.dt_left_obsA_nb = [idx_all.dt_left_obsA_nb; DT_30_left_obsA_NB_idx];
    idx_all.dt_right_obsA_nb = [idx_all.dt_right_obsA_nb; DT_30_right_obsA_NB_idx];
    
    idx_all.dt_left_obsB_bl = [idx_all.dt_left_obsB_bl; DT_30_left_obsB_block_idx];
    idx_all.dt_right_obsB_bl = [idx_all.dt_right_obsB_bl; DT_30_right_obsB_block_idx];
    idx_all.dt_left_obsB_nb = [idx_all.dt_left_obsB_nb; DT_30_left_obsB_NB_idx];
    idx_all.dt_right_obsB_nb = [idx_all.dt_right_obsB_nb; DT_30_right_obsB_NB_idx];
    
    idx_all.strt_obsA_left = [idx_all.strt_obsA_left; ST_30_strt_obsA_left_idx];
    idx_all.strt_obsB_left = [idx_all.strt_obsB_left; ST_30_strt_obsB_left_idx];
    idx_all.strt_obsA_right = [idx_all.strt_obsA_right; ST_30_strt_obsA_right_idx];
    idx_all.strt_obsB_right = [idx_all.strt_obsB_right; ST_30_strt_obsB_right_idx];
    
    idx_all.st_left_bln = [idx_all.st_left_bln; STT_30_left_baseline_idx];
    idx_all.st_right_bln = [idx_all.st_right_bln; STT_30_right_baseline_idx];
    idx_all.st_strt_bln = [idx_all.st_strt_bln; STT_strt_baseline_idx];
    
    idx_all.dt_bln = [idx_all.dt_bln;DTT_30_baseline_idx_bln] ;
    idx_all.dt_bln_left = [idx_all.dt_bln_left; DTT_30_bln_left_idx];
    idx_all.dt_bln_right = [idx_all.dt_bln_right; DTT_30_bln_right_idx];
    
    idx_all.dt_null_test = [idx_all.dt_null_test;DTT_30_test_idx];
    
    idx_all.dt_test_all = [idx_all.dt_test_all; DT_30_test_idx_all];
    idx_all.st_test_all = [idx_all.st_test_all; ST_30_test_idx_all];
    
    idx_all.strt_null = [idx_all.strt_null; STT_strt_null_idx];
    idx_all.strt_obs = [idx_all.strt_obs; STT_strt_obs_idx];
    
    idx_all.dt_all = [idx_all.dt_all; find(DT_idx_all(nn,:).*is30_all(nn,:))];
    st_block_all_idx = find(ST_idx_all(nn,:).*is30_all(nn,:).*is_obsA_all(nn,:).*is_obsBlock_all(nn,:));
    idx_all.st_all_obsA_bl{nn} = st_block_all_idx;
    st_nb_all_idx_tmp = find(ST_idx_all(nn,:).*is30_all(nn,:));
    st_nb_all_idx = st_nb_all_idx_tmp(~ismember(st_nb_all_idx_tmp,st_block_all_idx));
    idx_all.st_all_obsA_nb{nn} = st_nb_all_idx;
    
    %if nn==2, keyboard; end
    
    %% get the data
    
    %get STT null data
    xTraj_ST_30_left_null = xTraj_sub_all(nn,STT_30_left_baseline_idx);    yTraj_ST_30_left_null = yTraj_sub_all(nn,STT_30_left_baseline_idx);
    xTraj_ST_30_right_null = xTraj_sub_all(nn,STT_30_right_baseline_idx);    yTraj_ST_30_right_null = yTraj_sub_all(nn,STT_30_right_baseline_idx);
    xTraj_ST_strt_null = xTraj_sub_all(nn,STT_strt_baseline_idx);    yTraj_ST_strt_null = yTraj_sub_all(nn,STT_strt_baseline_idx);

    Time_ST_30_left_null=Time_sub_all(nn,STT_30_left_baseline_idx);
    Time_ST_30_right_null=Time_sub_all(nn,STT_30_right_baseline_idx);
    Time_ST_strt_null=Time_sub_all(nn,STT_strt_baseline_idx);
    
    %STT data, blocked
    xTraj_ST_30_leftA_block=xTraj_sub_all(nn,ST_30_left_obsA_block_idx);   yTraj_ST_30_leftA_block=yTraj_sub_all(nn,ST_30_left_obsA_block_idx);
    xTraj_ST_30_leftB_block=xTraj_sub_all(nn,ST_30_left_obsB_block_idx);   yTraj_ST_30_leftB_block=yTraj_sub_all(nn,ST_30_left_obsB_block_idx);
    
    xTraj_ST_30_rightA_block=xTraj_sub_all(nn,ST_30_right_obsA_block_idx);   yTraj_ST_30_rightA_block=yTraj_sub_all(nn,ST_30_right_obsA_block_idx);
    xTraj_ST_30_rightB_block=xTraj_sub_all(nn,ST_30_right_obsB_block_idx);   yTraj_ST_30_rightB_block=yTraj_sub_all(nn,ST_30_right_obsB_block_idx);
    
    Time_ST_30_leftA_block=Time_sub_all(nn,ST_30_left_obsA_block_idx);
    Time_ST_30_leftB_block=Time_sub_all(nn,ST_30_left_obsB_block_idx);
    
    Time_ST_30_rightA_block=Time_sub_all(nn,ST_30_right_obsA_block_idx);
    Time_ST_30_rightB_block=Time_sub_all(nn,ST_30_right_obsB_block_idx);
    
    %STT data, not blocked
    xTraj_ST_30_leftA_NB=xTraj_sub_all(nn,ST_30_left_obsA_NB_idx);   yTraj_ST_30_leftA_NB=yTraj_sub_all(nn,ST_30_left_obsA_NB_idx);
    xTraj_ST_30_leftB_NB=xTraj_sub_all(nn,ST_30_left_obsB_NB_idx);   yTraj_ST_30_leftB_NB=yTraj_sub_all(nn,ST_30_left_obsB_NB_idx);
    
    xTraj_ST_30_rightA_NB=xTraj_sub_all(nn,ST_30_right_obsA_NB_idx);   yTraj_ST_30_rightA_NB=yTraj_sub_all(nn,ST_30_right_obsA_NB_idx);
    xTraj_ST_30_rightB_NB=xTraj_sub_all(nn,ST_30_right_obsB_NB_idx);   yTraj_ST_30_rightB_NB=yTraj_sub_all(nn,ST_30_right_obsB_NB_idx);
    
    Time_ST_30_leftA_NB=Time_sub_all(nn,ST_30_left_obsA_NB_idx);
    Time_ST_30_leftB_NB=Time_sub_all(nn,ST_30_left_obsB_NB_idx);
    
    Time_ST_30_rightA_NB=Time_sub_all(nn,ST_30_right_obsA_NB_idx);
    Time_ST_30_rightB_NB=Time_sub_all(nn,ST_30_right_obsB_NB_idx);
    
    %baseline STT obstacle (focus only on blocked case for now)
    xTraj_ST_bln_30_leftA_block=xTraj_sub_all(nn,ST_bln_30_left_obsA_block_idx);   yTraj_ST_bln_30_leftA_block=yTraj_sub_all(nn,ST_bln_30_left_obsA_block_idx);
    xTraj_ST_bln_30_leftB_block=xTraj_sub_all(nn,ST_bln_30_left_obsB_block_idx);   yTraj_ST_bln_30_leftB_block=yTraj_sub_all(nn,ST_bln_30_left_obsB_block_idx);
    
    xTraj_ST_bln_30_rightA_block=xTraj_sub_all(nn,ST_bln_30_right_obsA_block_idx);   yTraj_ST_bln_30_rightA_block=yTraj_sub_all(nn,ST_bln_30_right_obsA_block_idx);
    xTraj_ST_bln_30_rightB_block=xTraj_sub_all(nn,ST_bln_30_right_obsB_block_idx);   yTraj_ST_bln_30_rightB_block=yTraj_sub_all(nn,ST_bln_30_right_obsB_block_idx);
    
    Time_ST_bln_30_leftA_block=Time_sub_all(nn,ST_bln_30_left_obsA_block_idx);
    Time_ST_bln_30_leftB_block=Time_sub_all(nn,ST_bln_30_left_obsB_block_idx);
    
    Time_ST_bln_30_rightA_block=Time_sub_all(nn,ST_bln_30_right_obsA_block_idx);
    Time_ST_bln_30_rightB_block=Time_sub_all(nn,ST_bln_30_right_obsB_block_idx);
    
    %Straight ahead case
    xTraj_strt_30_leftA = xTraj_sub_all(nn,ST_30_strt_obsA_left_idx);   yTraj_strt_30_leftA = yTraj_sub_all(nn,ST_30_strt_obsA_left_idx);
    xTraj_strt_30_leftB = xTraj_sub_all(nn,ST_30_strt_obsB_left_idx);   yTraj_strt_30_leftB = yTraj_sub_all(nn,ST_30_strt_obsB_left_idx);
    
    xTraj_strt_30_rightA = xTraj_sub_all(nn,ST_30_strt_obsA_right_idx);   yTraj_strt_30_rightA = yTraj_sub_all(nn,ST_30_strt_obsA_right_idx);
    xTraj_strt_30_rightB = xTraj_sub_all(nn,ST_30_strt_obsB_right_idx);   yTraj_strt_30_rightB = yTraj_sub_all(nn,ST_30_strt_obsB_right_idx);
    
    Time_strt_30_leftA=Time_sub_all(nn,ST_30_strt_obsA_left_idx);
    Time_strt_30_leftB=Time_sub_all(nn,ST_30_strt_obsB_left_idx);
    
    Time_strt_30_rightA=Time_sub_all(nn,ST_30_strt_obsA_right_idx);
    Time_strt_30_rightB=Time_sub_all(nn,ST_30_strt_obsB_right_idx);
    
    
    %DTT data
    xTraj_DT_30_left_blockA=xTraj_sub_all(nn,DT_30_left_obsA_block_idx);   yTraj_DT_30_left_blockA=yTraj_sub_all(nn,DT_30_left_obsA_block_idx);
    xTraj_DT_30_left_blockB=xTraj_sub_all(nn,DT_30_left_obsB_block_idx);   yTraj_DT_30_left_blockB=yTraj_sub_all(nn,DT_30_left_obsB_block_idx);
    
    xTraj_DT_30_right_blockA=xTraj_sub_all(nn,DT_30_right_obsA_block_idx);   yTraj_DT_30_right_blockA=yTraj_sub_all(nn,DT_30_right_obsA_block_idx);
    xTraj_DT_30_right_blockB=xTraj_sub_all(nn,DT_30_right_obsB_block_idx);   yTraj_DT_30_right_blockB=yTraj_sub_all(nn,DT_30_right_obsB_block_idx);
    
    Time_DT_30_left_blockA=Time_sub_all(nn,DT_30_left_obsA_block_idx);
    Time_DT_30_left_blockB=Time_sub_all(nn,DT_30_left_obsB_block_idx);
    
    Time_DT_30_right_blockA=Time_sub_all(nn,DT_30_right_obsA_block_idx);
    Time_DT_30_right_blockB=Time_sub_all(nn,DT_30_right_obsB_block_idx);
    
    %no block case
    xTraj_DT_30_left_NBA=xTraj_sub_all(nn,DT_30_left_obsA_NB_idx);   yTraj_DT_30_left_NBA=yTraj_sub_all(nn,DT_30_left_obsA_NB_idx);
    xTraj_DT_30_left_NBB=xTraj_sub_all(nn,DT_30_left_obsB_NB_idx);   yTraj_DT_30_left_NBB=yTraj_sub_all(nn,DT_30_left_obsB_NB_idx);
    
    xTraj_DT_30_right_NBA=xTraj_sub_all(nn,DT_30_right_obsA_NB_idx);   yTraj_DT_30_right_NBA=yTraj_sub_all(nn,DT_30_right_obsA_NB_idx);
    xTraj_DT_30_right_NBB=xTraj_sub_all(nn,DT_30_right_obsB_NB_idx);   yTraj_DT_30_right_NBB=yTraj_sub_all(nn,DT_30_right_obsB_NB_idx);
    
    
    Time_DT_30_left_NBA=Time_sub_all(nn,DT_30_left_obsA_NB_idx);
    Time_DT_30_left_NBB=Time_sub_all(nn,DT_30_left_obsB_NB_idx);
    
    Time_DT_30_right_NBA=Time_sub_all(nn,DT_30_right_obsA_NB_idx);
    Time_DT_30_right_NBB=Time_sub_all(nn,DT_30_right_obsB_NB_idx);
    
    %Lastly, interpolate the baseline DTT trials (want to plot this as well)
    xTraj_DT_30_left_bln = xTraj_sub_all(nn,DTT_30_bln_left_idx);   yTraj_DT_30_left_bln = yTraj_sub_all(nn,DTT_30_bln_left_idx);
    xTraj_DT_30_right_bln = xTraj_sub_all(nn,DTT_30_bln_right_idx);   yTraj_DT_30_right_bln = yTraj_sub_all(nn,DTT_30_bln_right_idx);
    
    Time_DT_30_left_bln = Time_sub_all(nn,DTT_30_bln_left_idx);
    Time_DT_30_right_bln = Time_sub_all(nn,DTT_30_bln_right_idx);
    
    
    %figure out which STT trials are optimal ones
    optimal_indicator_all(nn,:) = determine_path_type_821(xTraj_sub_all(nn,:),yTraj_sub_all(nn,:),ST_idx_all(nn,:),obs_idx_all(nn,:), ...
        isLeft_all(nn,:),isRight_all(nn,:), is30_all(nn,:),is30_all(nn,:) * nan );
    
    OO = optimal_indicator_all(nn,:);
    OO_ST_30_left_obsA = OO(ST_30_left_obsA_block_idx);
    OO_ST_30_right_obsA = OO(ST_30_right_obsA_block_idx);
    OO_ST_30_left_obsB = OO(ST_30_left_obsB_block_idx);
    OO_ST_30_right_obsB = OO(ST_30_right_obsB_block_idx);
    
    
    %% interpolate and window the data
    
    %ST 30, null
    [xTraj_ST_30_left_null,~]=window_data(Time_ST_30_left_null,xTraj_ST_30_left_null,trial_win);
    [yTraj_ST_30_left_null,~]=window_data(Time_ST_30_left_null,yTraj_ST_30_left_null,trial_win);
    
    [xTraj_ST_30_right_null,~]=window_data(Time_ST_30_right_null,xTraj_ST_30_right_null,trial_win);
    [yTraj_ST_30_right_null,~]=window_data(Time_ST_30_right_null,yTraj_ST_30_right_null,trial_win);
    
    [xTraj_ST_strt_null,~]=window_data(Time_ST_strt_null,xTraj_ST_strt_null,trial_win);
    [yTraj_ST_strt_null,~]=window_data(Time_ST_strt_null,yTraj_ST_strt_null,trial_win);
    
    %STT 30, blocked
    [xTraj_ST_30_leftA_block_binned,xTraj_ST_30_leftA_block_std]=window_data(Time_ST_30_leftA_block,xTraj_ST_30_leftA_block,trial_win);
    [xTraj_ST_30_leftB_block_binned,xTraj_ST_30_leftB_block_std]=window_data(Time_ST_30_leftB_block,xTraj_ST_30_leftB_block,trial_win);
    
    [yTraj_ST_30_leftA_block_binned,yTraj_ST_30_leftA_block_std]=window_data(Time_ST_30_leftA_block,yTraj_ST_30_leftA_block,trial_win);
    [yTraj_ST_30_leftB_block_binned,yTraj_ST_30_leftB_block_std]=window_data(Time_ST_30_leftB_block,yTraj_ST_30_leftB_block,trial_win);
    
    [xTraj_ST_30_rightA_block_binned,xTraj_ST_30_rightA_block_std]=window_data(Time_ST_30_rightA_block,xTraj_ST_30_rightA_block,trial_win);
    [xTraj_ST_30_rightB_block_binned,xTraj_ST_30_rightB_block_std]=window_data(Time_ST_30_rightB_block,xTraj_ST_30_rightB_block,trial_win);
    
    [yTraj_ST_30_rightA_block_binned,yTraj_ST_30_rightA_block_std]=window_data(Time_ST_30_rightA_block,yTraj_ST_30_rightA_block,trial_win);
    [yTraj_ST_30_rightB_block_binned,yTraj_ST_30_rightB_block_std]=window_data(Time_ST_30_rightB_block,yTraj_ST_30_rightB_block,trial_win);
    
    %STT 30, not blocked
    [xTraj_ST_30_leftA_NB_binned,xTraj_ST_30_leftA_NB_std]=window_data(Time_ST_30_leftA_NB,xTraj_ST_30_leftA_NB,trial_win);
    [xTraj_ST_30_leftB_NB_binned,xTraj_ST_30_leftB_NB_std]=window_data(Time_ST_30_leftB_NB,xTraj_ST_30_leftB_NB,trial_win);
    
    [yTraj_ST_30_leftA_NB_binned,yTraj_ST_30_leftA_NB_std]=window_data(Time_ST_30_leftA_NB,yTraj_ST_30_leftA_NB,trial_win);
    [yTraj_ST_30_leftB_NB_binned,yTraj_ST_30_leftB_NB_std]=window_data(Time_ST_30_leftB_NB,yTraj_ST_30_leftB_NB,trial_win);
    
    [xTraj_ST_30_rightA_NB_binned,xTraj_ST_30_rightA_NB_std]=window_data(Time_ST_30_rightA_NB,xTraj_ST_30_rightA_NB,trial_win);
    [xTraj_ST_30_rightB_NB_binned,xTraj_ST_30_rightB_NB_std]=window_data(Time_ST_30_rightB_NB,xTraj_ST_30_rightB_NB,trial_win);
    
    [yTraj_ST_30_rightA_NB_binned,yTraj_ST_30_rightA_NB_std]=window_data(Time_ST_30_rightA_NB,yTraj_ST_30_rightA_NB,trial_win);
    [yTraj_ST_30_rightB_NB_binned,yTraj_ST_30_rightB_NB_std]=window_data(Time_ST_30_rightB_NB,yTraj_ST_30_rightB_NB,trial_win);
    
    %STT w/ obs blocked (before test)
    [xTraj_ST_bln_30_leftA_block_binned,xTraj_ST_bln_30_leftA_block_std]=window_data(Time_ST_bln_30_leftA_block,xTraj_ST_bln_30_leftA_block,trial_win);
    [xTraj_ST_bln_30_leftB_block_binned,xTraj_ST_bln_30_leftB_block_std]=window_data(Time_ST_bln_30_leftB_block,xTraj_ST_bln_30_leftB_block,trial_win);
    
    [yTraj_ST_bln_30_leftA_block_binned,yTraj_ST_bln_30_leftA_block_std]=window_data(Time_ST_bln_30_leftA_block,yTraj_ST_bln_30_leftA_block,trial_win);
    [yTraj_ST_bln_30_leftB_block_binned,yTraj_ST_bln_30_leftB_block_std]=window_data(Time_ST_bln_30_leftB_block,yTraj_ST_bln_30_leftB_block,trial_win);
    
    [xTraj_ST_bln_30_rightA_block_binned,xTraj_ST_bln_30_rightA_block_std]=window_data(Time_ST_bln_30_rightA_block,xTraj_ST_bln_30_rightA_block,trial_win);
    [xTraj_ST_bln_30_rightB_block_binned,xTraj_ST_bln_30_rightB_block_std]=window_data(Time_ST_bln_30_rightB_block,xTraj_ST_bln_30_rightB_block,trial_win);
    
    [yTraj_ST_bln_30_rightA_block_binned,yTraj_ST_bln_30_rightA_block_std]=window_data(Time_ST_bln_30_rightA_block,yTraj_ST_bln_30_rightA_block,trial_win);
    [yTraj_ST_bln_30_rightB_block_binned,yTraj_ST_bln_30_rightB_block_std]=window_data(Time_ST_bln_30_rightB_block,yTraj_ST_bln_30_rightB_block,trial_win);
    
    
    %Straight 30
    [xTraj_strt_30_leftA_binned,xTraj_strt_30_leftA_std]=window_data(Time_strt_30_leftA,xTraj_strt_30_leftA,trial_win);
    [yTraj_strt_30_leftA_binned,yTraj_strt_30_leftA_std]=window_data(Time_strt_30_leftA,yTraj_strt_30_leftA,trial_win);
    
    [xTraj_strt_30_rightA_binned,xTraj_strt_30_rightA_std]=window_data(Time_strt_30_rightA,xTraj_strt_30_rightA,trial_win);
    [yTraj_strt_30_rightA_binned,yTraj_strt_30_rightA_std]=window_data(Time_strt_30_rightA,yTraj_strt_30_rightA,trial_win);
    
    [xTraj_strt_30_leftB_binned,xTraj_strt_30_leftB_std]=window_data(Time_strt_30_leftB,xTraj_strt_30_leftB,trial_win);
    [yTraj_strt_30_leftB_binned,yTraj_strt_30_leftB_std]=window_data(Time_strt_30_leftB,yTraj_strt_30_leftB,trial_win);
    
    [xTraj_strt_30_rightB_binned,xTraj_strt_30_rightB_std]=window_data(Time_strt_30_rightB,xTraj_strt_30_rightB,trial_win);
    [yTraj_strt_30_rightB_binned,yTraj_strt_30_rightB_std]=window_data(Time_strt_30_rightB,yTraj_strt_30_rightB,trial_win);
    
    
    %DTT 30 blocked
    [xTraj_DT_30_left_blockA_binned,xTraj_ST_30_left_blockA_std]=window_data(Time_DT_30_left_blockA,xTraj_DT_30_left_blockA,trial_win);
    [xTraj_DT_30_left_blockB_binned,xTraj_ST_30_left_blockB_std]=window_data(Time_DT_30_left_blockB,xTraj_DT_30_left_blockB,trial_win);
    
    [yTraj_DT_30_left_blockA_binned,yTraj_ST_30_left_blockA_std]=window_data(Time_DT_30_left_blockA,yTraj_DT_30_left_blockA,trial_win);
    [yTraj_DT_30_left_blockB_binned,yTraj_ST_30_left_blockB_std]=window_data(Time_DT_30_left_blockB,yTraj_DT_30_left_blockB,trial_win);
    
    [xTraj_DT_30_right_blockA_binned,xTraj_ST_30_right_blockA_std]=window_data(Time_DT_30_right_blockA,xTraj_DT_30_right_blockA,trial_win);
    [xTraj_DT_30_right_blockB_binned,xTraj_ST_30_right_blockB_std]=window_data(Time_DT_30_right_blockB,xTraj_DT_30_right_blockB,trial_win);
    
    [yTraj_DT_30_right_blockA_binned,yTraj_ST_30_right_blockA_std]=window_data(Time_DT_30_right_blockA,yTraj_DT_30_right_blockA,trial_win);
    [yTraj_DT_30_right_blockB_binned,yTraj_ST_30_right_blockB_std]=window_data(Time_DT_30_right_blockB,yTraj_DT_30_right_blockB,trial_win);
    
    %DTT 30 NB (not blocked)
    [xTraj_DT_30_left_NBA_binned,xTraj_ST_30_left_NBA_std]=window_data(Time_DT_30_left_NBA,xTraj_DT_30_left_NBA,trial_win);
    [xTraj_DT_30_left_NBB_binned,xTraj_ST_30_left_NBB_std]=window_data(Time_DT_30_left_NBB,xTraj_DT_30_left_NBB,trial_win);
    
    [yTraj_DT_30_left_NBA_binned,yTraj_ST_30_left_NBA_std]=window_data(Time_DT_30_left_NBA,yTraj_DT_30_left_NBA,trial_win);
    [yTraj_DT_30_left_NBB_binned,yTraj_ST_30_left_NBB_std]=window_data(Time_DT_30_left_NBB,yTraj_DT_30_left_NBB,trial_win);
    
    [xTraj_DT_30_right_NBA_binned,xTraj_ST_30_right_NBA_std]=window_data(Time_DT_30_right_NBA,xTraj_DT_30_right_NBA,trial_win);
    [xTraj_DT_30_right_NBB_binned,xTraj_ST_30_right_NBB_std]=window_data(Time_DT_30_right_NBB,xTraj_DT_30_right_NBB,trial_win);
    
    [yTraj_DT_30_right_NBA_binned,yTraj_ST_30_right_NBA_std]=window_data(Time_DT_30_right_NBA,yTraj_DT_30_right_NBA,trial_win);
    [yTraj_DT_30_right_NBB_binned,yTraj_ST_30_right_NBB_std]=window_data(Time_DT_30_right_NBB,yTraj_DT_30_right_NBB,trial_win);
    
    
    %DT 30 baseline (left and right)
    [xTraj_DT_30_left_bln_binned,~] = window_data(Time_DT_30_left_bln, xTraj_DT_30_left_bln, trial_win);
    [yTraj_DT_30_left_bln_binned,~] = window_data(Time_DT_30_left_bln, yTraj_DT_30_left_bln, trial_win);
    
    [xTraj_DT_30_right_bln_binned,~] = window_data(Time_DT_30_right_bln, xTraj_DT_30_right_bln, trial_win);
    [yTraj_DT_30_right_bln_binned,~] = window_data(Time_DT_30_right_bln, yTraj_DT_30_right_bln, trial_win);
    
    %%%%% Lets combine the blocked and no blocked cases (left and right separately)
%     Time_DT_30_obs_leftA = [Time_DT_30_left_blockA, Time_DT_30_right_NBA];
%     xTraj_DT_30_obs_leftA = [xTraj_DT_30_left_blockA,xTraj_DT_30_right_NBA ];
%     yTraj_DT_30_obs_leftA = [yTraj_DT_30_left_blockA,yTraj_DT_30_right_NBA];
%     
%     Time_DT_30_obs_rightA = [Time_DT_30_right_blockA, Time_DT_30_left_NBA];
%     xTraj_DT_30_obs_rightA = [xTraj_DT_30_right_blockA,xTraj_DT_30_left_NBA ];
%     yTraj_DT_30_obs_rightA = [yTraj_DT_30_right_blockA,yTraj_DT_30_left_NBA];
    
    %now interpolate
%     [xTraj_DT_30_obs_leftA_binned,~]=window_data(Time_DT_30_obs_leftA,xTraj_DT_30_obs_leftA,trial_win);
%     [yTraj_DT_30_obs_leftA_binned,~]=window_data(Time_DT_30_obs_leftA,yTraj_DT_30_obs_leftA,trial_win);
%     
%     [xTraj_DT_30_obs_rightA_binned,~]=window_data(Time_DT_30_obs_rightA,xTraj_DT_30_obs_rightA,trial_win);
%     [yTraj_DT_30_obs_rightA_binned,~]=window_data(Time_DT_30_obs_rightA,yTraj_DT_30_obs_rightA,trial_win);
    
    
    %% get the success rate
    
    %IMPORTANT: might want to get rid of the glass criteria in bad labels here
    %     sub_good_trial = ~sub_bad_trial;
    %     MT_label_sub = MT_labels_sub_all(nn,:)'; %note that bad label doesnt include MT
    %     MT = MT_sub_all(nn,:);
    %
    %     good_label= prod(sub_good_trial,2) .*MT_label_sub; %now we have overall good trials
    
    sub_rwd_label = rwd_label(nn,:);
    
    
    %now get it for all cases
    
    %get it across target directions for now
    good_st_obsA_bl =sub_rwd_label( [ST_30_left_obsA_block_idx, ST_30_right_obsA_block_idx]);
    good_st_obsA_nb = sub_rwd_label( [ST_30_left_obsA_NB_idx, ST_30_right_obsA_NB_idx]);
    good_st_obsB_bl =sub_rwd_label( [ST_30_left_obsB_block_idx, ST_30_right_obsB_block_idx]);
    good_st_obsB_nb = sub_rwd_label( [ST_30_left_obsB_NB_idx, ST_30_right_obsB_NB_idx]);
    
    good_dt_obsA_bl =sub_rwd_label( [DT_30_left_obsA_block_idx, DT_30_right_obsA_block_idx]);
    good_dt_obsA_nb = sub_rwd_label( [DT_30_left_obsA_NB_idx, DT_30_right_obsA_NB_idx]);
    good_dt_obsB_bl =sub_rwd_label( [DT_30_left_obsB_block_idx, DT_30_right_obsB_block_idx]);
    good_dt_obsB_nb = sub_rwd_label( [DT_30_left_obsB_NB_idx, DT_30_right_obsB_NB_idx]);
    
    
    %calculate the percentages
    score_st_obsA_bl = sum(good_st_obsA_bl) / length( good_st_obsA_bl ) * 100;
    score_st_obsA_nb = sum(good_st_obsA_nb) / length( good_st_obsA_nb ) * 100;
    score_st_obsB_bl = sum(good_st_obsB_bl) / length( good_st_obsB_bl ) * 100;
    score_st_obsB_nb = sum(good_st_obsB_nb) / length( good_st_obsB_nb ) * 100;
    
    score_dt_obsA_bl = sum(good_dt_obsA_bl) / length( good_dt_obsA_bl ) * 100;
    score_dt_obsA_nb = sum(good_dt_obsA_nb) / length( good_dt_obsA_nb ) * 100;
    score_dt_obsB_bl = sum(good_dt_obsB_bl) / length( good_dt_obsB_bl ) * 100;
    score_dt_obsB_nb = sum(good_dt_obsB_nb) / length( good_dt_obsB_nb ) * 100;
    
    %% plotting
    
    %IMPORTANT: See previous versions for subplots of all the data
    %if nn==num_subjects, close all; end
    
    if show_fig_flag
        
        h1=figure(fig_num); %STT blocked data during test
        
        %%%Update (5/30): Instead of superimposing the left and right cases
        %%%on the same plot, separate them, and then also show the movement
        %%%to the opposite target as well
        
%         %ST 30 A
%         subplot(ceil(num_subjects/4),4,nn); %lets put each subject on an individual subplot
%         
%         display_obs_plot_V2(xTraj_ST_30_leftA_block_binned,yTraj_ST_30_leftA_block_binned,linewidth,light_blue,OO_ST_30_left_obsA);
%         display_obs_plot_V2(xTraj_ST_30_rightA_block_binned,yTraj_ST_30_rightA_block_binned,linewidth,red_,OO_ST_30_right_obsA);
%         
%         plot(startx,starty,'color',start_col,'linewidth',2.5)
%         plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',2.5,'linestyle','--','color',light_blue)
%         plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',2.5,'linestyle','--','color',red_)
%         plot(LT30_obs{1}(:,1),LT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
%         plot(RT30_obs{1}(:,1),RT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
%         patch(LT30_obs{1}(:,1),LT30_obs{1}(:,2),light_blue)
%         patch(RT30_obs{1}(:,1),RT30_obs{1}(:,2),red_)
%         
%         %title(['Score: ', num2str(score_st_obsA_bl), '%']);
%         axis equal;
%         xlabel(Sub_id_all(nn),'FontWeight','bold');
%         set(gca,'ytick','','xtick','');
%         %axis([5 40 25 50])
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %ST 30 A
        subplot(ceil(num_subjects/4),4,nn); %lets put each subject on an individual subplot
        hold on;
        
        display_obs_plot_V2(xTraj_ST_30_rightA_block_binned,yTraj_ST_30_rightA_block_binned,linewidth,red_,OO_ST_30_right_obsA);
        display_obs_plot_V2(xTraj_ST_30_leftA_NB_binned, yTraj_ST_30_leftA_NB_binned,linewidth, light_blue, ...
            ones(size(xTraj_ST_30_leftA_NB_binned,2),1) );
        
        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot(RT30_obs{1}(:,1),RT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        patch(RT30_obs{1}(:,1),RT30_obs{1}(:,2),grey)
        
        %title(['Score: ', num2str(score_st_obsA_bl), '%']);
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        %axis([5 40 25 50])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        h2=figure(fig_num+1);
        
        %DT 30 A (for now just plot right obstacle data)
        subplot(ceil(num_subjects/4),4,nn); hold on;
        %xTraj_DT_30_left_blockA_binned
        display_obs_plot(xTraj_DT_30_right_blockA_binned,yTraj_DT_30_right_blockA_binned,linewidth,red_,linestyle1); %these are combined data
        display_obs_plot(xTraj_DT_30_left_NBA_binned,yTraj_DT_30_left_NBA_binned,linewidth,light_blue,linestyle1);
        
        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        %plot(LT30_obs{1}(:,1),LT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        plot(RT30_obs{1}(:,1),RT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        %patch(LT30_obs{1}(:,1),LT30_obs{1}(:,2),light_blue)
        patch(RT30_obs{1}(:,1),RT30_obs{1}(:,2),grey)
        
        %title(['Score: ', num2str(score_dt_obsA_bl), '%']);
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %straight ahead cases, lets do 2 different figures for each obstacle (Left & right)
        h3 = figure(fig_num+2);
        
        subplot(ceil(num_subjects/4),4,nn); hold on;
        
        display_obs_plot(xTraj_strt_30_leftA_binned,yTraj_strt_30_leftA_binned,linewidth,black_,linestyle1);
        
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k');
        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot(LT30_obs{1}(:,1),LT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        patch(LT30_obs{1}(:,1),LT30_obs{1}(:,2),grey)
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        %if rr==1,title('Obs A, Left, 30 Deg');end;        
        
        
        h4 = figure(fig_num+3);
        subplot(ceil(num_subjects/4),4,nn); hold on;
        
        display_obs_plot(xTraj_strt_30_rightA_binned,yTraj_strt_30_rightA_binned,linewidth,black_,linestyle1);
        
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot(RT30_obs{1}(:,1),RT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        patch(RT30_obs{1}(:,1),RT30_obs{1}(:,2),grey)
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        %if rr==2, title('Obs A, Right, 30 Deg');end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Also plot null 2-target cases
        h5 = figure(fig_num+4);
        subplot(ceil(num_subjects/4),4,nn); hold on;
        
        display_obs_plot(xTraj_DT_30_left_bln_binned,yTraj_DT_30_left_bln_binned,linewidth,light_blue,linestyle1);
        display_obs_plot(xTraj_DT_30_right_bln_binned,yTraj_DT_30_right_bln_binned,linewidth,red_,linestyle1);

        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %get 2 tgt left obstacle data as well
        h6=figure(fig_num+5);
        
        subplot(ceil(num_subjects/4),4,nn);  hold on;

        display_obs_plot(xTraj_DT_30_right_NBA_binned,yTraj_DT_30_right_NBA_binned,linewidth,red_,linestyle1);
        display_obs_plot(xTraj_DT_30_left_blockA_binned,yTraj_DT_30_left_blockA_binned,linewidth,light_blue,linestyle1);
        
        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot(LT30_obs{1}(:,1),LT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        %plot(RT30_obs{1}(:,1),RT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        patch(LT30_obs{1}(:,1),LT30_obs{1}(:,2),grey)
        %patch(RT30_obs{1}(:,1),RT30_obs{1}(:,2),red_)
        
        %title(['Score: ', num2str(score_dt_obsA_bl), '%']);
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        
        
        %plot STT with left blocked data
        h7 = figure(fig_num+6);
        subplot(ceil(num_subjects/4),4,nn); hold on;
        
        display_obs_plot_V2(xTraj_ST_30_leftA_block_binned,yTraj_ST_30_leftA_block_binned,linewidth,light_blue,...
            OO_ST_30_left_obsA);
        display_obs_plot_V2(xTraj_ST_30_rightA_NB_binned, yTraj_ST_30_rightA_NB_binned,linewidth, red_,...
            ones(size(xTraj_ST_30_rightA_NB_binned,1)) );
        
        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot(LT30_obs{1}(:,1),LT30_obs{1}(:,2),'color',grey,'linewidth',1.5)
        patch(LT30_obs{1}(:,1),LT30_obs{1}(:,2),grey)
        
        %title(['Score: ', num2str(score_st_obsA_bl), '%']);
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        %axis([5 40 25 50])
        
        %plot STT null left and right data
        h8 = figure(fig_num+7);
        subplot(ceil(num_subjects/4),4,nn); hold on;
        
        plot(xTraj_ST_30_left_null,yTraj_ST_30_left_null,'linewidth', linewidth,'color', 'b');
        plot(xTraj_ST_30_right_null,yTraj_ST_30_right_null,'linewidth', linewidth,'color', 'r');
        plot(xTraj_ST_strt_null,yTraj_ST_strt_null,'linewidth', linewidth,'color', grey);

        plot(startx,starty,'color',start_col,'linewidth',2.5)
        plot([tsx,left_30_tgt(1)],[tsy,left_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        plot([tsx,right_30_tgt(1)],[tsy,right_30_tgt(2)],'linewidth',1,'linestyle','--','color','k')
        
        axis equal;
        xlabel(Sub_id_all(nn),'FontWeight','bold');
        set(gca,'ytick','','xtick','');
        
    end
    
    
end

%% add titles

% if show_fig_flag
%     
%     H=[h1,h2,h3,h4,h5,h6,h7];
%     %obs_text={['Obstacle B'], ['Obstacle A']};
%     titles={['Single Target Trials, 30 Degrees w/ obsR (blocked)'], ['Double Target Trials, w/ obsR, 30 Degrees'], ['Straight ahead, left obs'],...
%         ['Straight ahead, right obs'], ['Double Target Trials, null'], ['Double Target Trials, w/ obsL, 30 Degrees'],['Single Target Trials, 30 Degrees w/ obsL (blocked)']};
%     for tt=1:length(H)
%         figure(fig_num+(tt-1));
%         axes_h=findobj(H(tt), 'Type', 'Axes');
%         ylabel_h=get(axes_h, 'YLabel');
%         xlabel_h=get(axes_h, 'XLabel');
%         overall_ylabel=[ylabel_h{2:2:end}];
%         overall_xlabel=[xlabel_h{1:2}];
%         
% %         if tt~= length(H)
% %             for qq=1:length(overall_ylabel)
% %                 set(overall_ylabel(qq), 'String', Sub_id_all((end+1)-(qq)),'FontWeight','bold') %double check to make sure sub id's are corect on plot
% %             end
% %             
% %             for qq=1:length(overall_xlabel)
% %                 set(overall_xlabel(qq), 'String', obs_text(qq),'FontWeight','bold')
% %             end
% %         end
%         
%         suptitle(titles(tt));
%         
%         set(axes_h, 'Fontsize', 14);
%         %text(-20,125,titles(tt), 'FontSize',16);
%     end
%     
% end


%%
%keyboard;
%save our variables to use in other analyses
save('group_analysis.mat');













































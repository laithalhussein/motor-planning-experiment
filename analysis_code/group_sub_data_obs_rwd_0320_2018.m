function [dat_all,info_all]=group_sub_data_obs_rwd_0320_2018(Cdr)
%remember that we only have 2 types of obstacles for the latest version of the experiment
close all;
clearvars -except Cdr
home

dat_filename=ls('obs_exp_parse_all*'); 
all_data=load(strcat(Cdr,'\',dat_filename));
all_data=all_data.DAT; 
num_subjects=length(all_data);

%Exp sequence (forward direction only) (should be 15 blocks total)

%fam (60 trials) -> fam (60 trials) -> ST obs (60 trials) -> ST obs (60 trials) -> DTT + STT baseline (60 trials) x 2
%-> ST obs (60 trials) -> DTT + STT baseline (60 trials) -> ST obs (60 trials) ->
%obs 2 cm test (100 trials) * 4 -> obs 1 cm test (100 trials) * 2

trials_per_block=[(cellfun(@(x) length(x),all_data{1}.tgt_file)/2)-1]';
Nb = length(trials_per_block); %total number of blocks

%New experiment sequence:
% 2 Fam (60)-> baseline stt w/ obs x 2 (50) -> baseline dtt + stt x 2 (60) -> baseline stt w/ obs (30) ->
% baseline dtt+ stt (60) -> baseline stt w/ obs (30) -> TEST!!!
%TEST: A-A-A x 2...B-B-B x 2 (2 for left and right obstacles)

exp_seq={ [1:120],[121:220],[221:340],[341:370],[371:430], [431:460], [461:760], [761:1060] };

%%%%%%%%%%%Preallocation%%%%%%%%
%Main Data - also add ting time(?)
total_forward_trials=sum(trials_per_block);
xTraj_sub_all=cell(num_subjects,total_forward_trials);
yTraj_sub_all=cell(num_subjects,total_forward_trials);
xTraj_screen_sub_all=cell(num_subjects,total_forward_trials);
yTraj_screen_sub_all=cell(num_subjects,total_forward_trials);
Vel_trav_sub_all=cell(num_subjects,total_forward_trials);
Pos_sub_all=cell(num_subjects,total_forward_trials);
Time_sub_all=cell(num_subjects,total_forward_trials);
Time_vec_sub_all=cell(num_subjects,total_forward_trials);
MT_sub_all=zeros(num_subjects,total_forward_trials);
react_time_mo_sub_all=zeros(num_subjects,total_forward_trials);
Score_sub_all=zeros(num_subjects,length(all_data{1}.tgt_file));
bad_labels_sub_all=zeros(num_subjects,total_forward_trials,3);
MT_labels_sub_all=zeros(num_subjects,total_forward_trials);
rwd_fb_sub_all = zeros(num_subjects,total_forward_trials);
Sub_id_all=cell(num_subjects,1);
exp_type_all = nan(num_subjects,1);

Ideal_dir_all=zeros(num_subjects,total_forward_trials);
ST_idx_all=zeros(num_subjects,total_forward_trials);   DT_idx_all=zeros(num_subjects,total_forward_trials);
isLeft_all=zeros(num_subjects,total_forward_trials);   isRight_all=zeros(num_subjects,total_forward_trials);   isStraight_all=zeros(num_subjects,total_forward_trials);
obs_idx_all=zeros(num_subjects,total_forward_trials);
is_obsA_all=zeros(num_subjects,total_forward_trials);   is_obsB_all=zeros(num_subjects,total_forward_trials); 
is30_all=zeros(num_subjects,total_forward_trials);   is45_all=zeros(num_subjects,total_forward_trials);   is90_all=zeros(num_subjects,total_forward_trials);
is_obsLeft_all=zeros(num_subjects,total_forward_trials);   is_obsRight_all=zeros(num_subjects,total_forward_trials);
is_obsBlock_all=zeros(num_subjects,total_forward_trials);   is_obsNB_all=zeros(num_subjects,total_forward_trials);

index_per_block = cell(Nb,1);

for n=1:num_subjects
    
    dat=all_data{n};
    traj=dat.MovTab;
    traj_screen=dat.Mov;
    vel_trav=dat.Vel_trav;
    time=dat.Time_tab;
    time_vec=dat.Time_vec;
    pos=dat.Displ;
    mt=dat.MT;
    bad_labels=vertcat(dat.bad_labels{:});
    mt_labels=vertcat(dat.MT_labels{:});
    react_time_mo=dat.react_time_mo;
    score=dat.Score;
    sub_id=dat.sub_id;
    rwd_fb = dat.rwd_fb;
    
    %tgt files
    tgt_file=dat.tgt_file;
    tgt_all=vertcat(tgt_file{:});
    
    %intended direction of movements -DONT use this for determining right vs left targets
    ideal_dir=dat.Ideal_dir; %Somthing is wrong is Ideal_dir. The angle might be calculated correctly but the direction (sign) is not.
    
    %extract index where data was actually collected
    start_idx = 2;
%     forward_start_idx = 2;
    for kk=1:Nb
        
        current_idx = [start_idx : length(vertcat(tgt_file{1:kk})) - 1];
        index_per_block{kk} = current_idx;
        start_idx = current_idx(end) + 3;
        
%         forward_idx_tmp = current_idx(1:2:end);
%         if kk> 1, forward_start_idx = forward_idx_tmp(end) + 4; end
%         if tgt_flag==1, exp_seq{kk} = [forward_start_idx: 2: length(vertcat(tgt_file{1:kk})) - 1]; end
    end
    
    relevant_trials_idx=[index_per_block{:}];
    relevant_tgt=tgt_all(relevant_trials_idx,:); %get the relevant tgt files
    
    forward_trials_idx=find(relevant_tgt(:,4)~=975); %975 is the only unique coordinate for the starting target
    forward_tgt=relevant_tgt(forward_trials_idx,:);
    ideal_dir_forward=ideal_dir(abs(ideal_dir)<= 90);
    mov_angles=unique(abs(ideal_dir_forward),'sorted');
    
    %single and double target trials
    ST_idx=(forward_tgt(:,11)==1)';
    DT_idx=(forward_tgt(:,11)==3)';
    
    %is the target left, right or straight ahead
    isLeft=zeros(1,length(ideal_dir_forward)); %left, right, or straight ahead target
    isRight=isLeft; isStraight=isLeft;
    for kk=1:length(forward_tgt)
        if forward_tgt(kk,11)==1
            if forward_tgt(kk,3) > 950
                isRight(kk)=1;
            elseif forward_tgt(kk,3) < 950
                isLeft(kk)=1;
            elseif forward_tgt(kk,3) == 950
                isStraight(kk)=1;
            end
        elseif forward_tgt(kk,11)==3
            if forward_tgt(kk,12)==0
                isLeft(kk)=1;
            elseif forward_tgt(kk,12)==1
                isRight(kk)=1;
            end
        end
    end
    %isRight(ideal_dir_forward<0)=1; isLeft(ideal_dir_forward>0)=1; isStraight(ideal_dir_forward==0)=1;
    
    obs_tgt_idx=[forward_tgt(:,13)>0]';%is there an obstacle (1) or not (0)
    %different obstacle types -- in this case let A be 2cm and B be 1 cm
    is_obsA=zeros(1,length(obs_tgt_idx));
    is_obsB=is_obsA;
    
    is_obsA(forward_tgt(:,13)==2 & forward_tgt(:,15)==2)=1; 
    is_obsB(forward_tgt(:,13)==2 & forward_tgt(:,15)==1)=1;
    
    %angle of the movement
    is30=0*ideal_dir_forward;
    is90=is30;
    
    is90(abs(ideal_dir_forward)==mov_angles(1))=1;      is30(abs(ideal_dir_forward)==mov_angles(2))=1; 
    
    %is the obstacle on the left or right
    obs_pos=[0*forward_tgt(:,14)]';
    is_obsLeft=obs_pos; is_obsRight=obs_pos;
    is_obsLeft(forward_tgt(:,13)>0 & forward_tgt(:,14)==0)=1;
    is_obsRight(forward_tgt(:,13)>0 & forward_tgt(:,14)==1)=1;
    
    %is the obstacle blocking the target being cued-check for STT *and* DTT
    obs_block=[0*forward_tgt(:,14)]';
    is_obsBlock=obs_block; is_obsNB=obs_block;
    
%     is_obsBlock(forward_tgt(:,11)==3 & forward_tgt(:,13)>0 & forward_tgt(:,14)==forward_tgt(:,12))=1;
%     is_obsNB(forward_tgt(:,11)==3 & forward_tgt(:,13)>0 & forward_tgt(:,14)~=forward_tgt(:,12))=1;

    is_obsBlock( (obs_tgt_idx & (isLeft == is_obsLeft)) | (obs_tgt_idx & isRight == is_obsRight) ) = 1;
    is_obsNB( (obs_tgt_idx & isLeft ~= is_obsLeft) | (obs_tgt_idx & isRight ~= is_obsRight) ) = 1;
    
    %start saving data
    Ideal_dir_all(n,:)=ideal_dir_forward;
    ST_idx_all(n,:)=ST_idx;   DT_idx_all(n,:)=DT_idx;
    isLeft_all(n,:)=isLeft;   isRight_all(n,:)=isRight;   isStraight_all(n,:)=isStraight;
    obs_idx_all(n,:)=obs_tgt_idx;
    is_obsA_all(n,:)=is_obsA;   is_obsB_all(n,:)=is_obsB; 
    is30_all(n,:)=is30;    is90_all(n,:)=is90;
    is_obsLeft_all(n,:)=is_obsLeft;   is_obsRight_all(n,:)=is_obsRight;
    is_obsBlock_all(n,:)=is_obsBlock;
    is_obsNB_all(n,:)=is_obsNB;
    
    %IMPORTANT: Filtering should happen here

    for k=1:length(forward_trials_idx) 
        trial_idx=forward_trials_idx(k);
        xTraj_sub_all{n,k}=traj{trial_idx}(:,1);
        yTraj_sub_all{n,k}=traj{trial_idx}(:,2);
        xTraj_screen_sub_all{n,k}=traj_screen{trial_idx}(:,1);
        yTraj_screen_sub_all{n,k}=traj_screen{trial_idx}(:,2);
        Vel_trav_sub_all{n,k}=vel_trav{trial_idx};
        Pos_sub_all{n,k}=pos{trial_idx};
        Time_sub_all{n,k}=time{trial_idx};
        Time_vec_sub_all{n,k}=time_vec{trial_idx};
    end
    
    MT_sub_all(n,:)=mt(forward_trials_idx);
    react_time_mo_sub_all(n,:)=react_time_mo(forward_trials_idx);
    bad_labels_sub_all(n,:,:)=bad_labels(forward_trials_idx,:);
    MT_labels_sub_all(n,:)=mt_labels(forward_trials_idx);
    rwd_fb_sub_all(n,:) = rwd_fb(forward_trials_idx);
    
    norm_score=(score./trials_per_block)*100;
    Score_sub_all(n,:)=norm_score;
    Sub_id_all{n}=sub_id;
    
    %save which experiment version was tested for each subject
    first_obsA_idx = find(is_obsA==1,1,'first');
    first_obsB_idx = find(is_obsB==1,1,'first');
    
    if first_obsA_idx < first_obsB_idx
        exp_type_all(n) = 1;
    else
        exp_type_all(n) = 0;
    end
    
    %keyboard;
    
end

%saving the data
info_all.Sub_id_all=Sub_id_all;
info_all.Score_sub_all=Score_sub_all;
info_all.MT_labels_sub_all=MT_labels_sub_all;
info_all.Bad_labels_sub_all=bad_labels_sub_all;
info_all.forward_tgt = forward_tgt;
info_all.num_subjects=num_subjects;
info_all.Ideal_dir_all=Ideal_dir_all;
info_all.exp_seq=exp_seq;   info_all.trials_per_block=trials_per_block;
info_all.ST_idx_all=ST_idx_all;   info_all.DT_idx_all=DT_idx_all;
info_all.isLeft_all=isLeft_all;   info_all.isRight_all=isRight_all;   info_all.isStraight_all=isStraight_all;
info_all.obs_idx_all=obs_idx_all;
info_all.is_obsA_all=is_obsA_all;   info_all.is_obsB_all=is_obsB_all; 
info_all.is30_all=is30_all;   info_all.is90_all=is90_all;
info_all.is_obsLeft_all=is_obsLeft_all;   info_all.is_obsRight_all=is_obsRight_all;
info_all.is_obsBlock_all=is_obsBlock_all;   info_all.is_obsNB_all=is_obsNB_all;
info_all.tgt = forward_tgt;
info_all.rwd_fb_sub_all = rwd_fb_sub_all;
info_all.exp_type = exp_type_all;

dat_all.xTraj_sub_all=xTraj_sub_all;
dat_all.yTraj_sub_all=yTraj_sub_all;
dat_all.xTraj_screen_sub_all=xTraj_screen_sub_all;
dat_all.yTraj_screen_sub_all=yTraj_screen_sub_all;
dat_all.Pos_sub_all=Pos_sub_all;
dat_all.Vel_trav_sub_all=Vel_trav_sub_all;
dat_all.react_time_mo_sub_all=react_time_mo_sub_all;
dat_all.Time_sub_all=Time_sub_all;
dat_all.Time_vec_sub_all=Time_vec_sub_all;
dat_all.MT_sub_all=MT_sub_all;

save(['obs_exp_dat_all_',num2str(num_subjects),'_subjects'],'dat_all','info_all');

end




    
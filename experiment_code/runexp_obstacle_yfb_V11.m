function runexp_obstacle_yfb_V11(name_prefix,tgt_file_name_prefix,tgt_set)
close all;

%%NOTE: could probably parallelize a lot of this....

%% Initialization 
% Main screen is on the right, secondary(tester) screen on the left
screenNumber = 0;
Screen('Preference','EmulateOldPTB',0);
Screen('Preference','TextRenderer',0);

%%% Load sound %%%
load ting
load glass %fs1
load obstacle_error
load mario_coin %fs2
load car_alarm %fs5
load MT_bad_buzzer %fs6
ting=ting'; %#ok<NODEF>
ting=ting(:,1:17640);
glass=glass'; %#ok<NODEF>
obstacle_error=obstacle_error';
obstacle_error=obstacle_error(:,3400:50000);
mario_coin=mario_coin';
mario_coin=mario_coin(:,15000:50000);
car_alarm=car_alarm';
car_alarm=car_alarm(:,(3400:3400+20000));
MT_bad_buzzer=MT_bad_buzzer';
MT_bad_buzzer=MT_bad_buzzer(:,30000:30000+10000);

desired_ting_time=0.1;
ting_samples=desired_ting_time*fs;
ting=ting(:,6000:6000+ting_samples); %tmake the ting last the desired ting time

ting_time=(1/fs)*length(ting);

grey=[256/2,256/2,256/2];

% tgt_file_name_prefix='obs_exp_45';
% name_prefix='tst';
% tgt_set='a';

%%% Load target file %%%
psychtest=load([tgt_file_name_prefix,tgt_set,'.tgt']); 

% Screen is able to do a lot of configuration and performance checks on
% open, and will print out a fair amount of detailed information when
% it does.  These commands supress that checking behavior and just let
% the demo go straight into action.  See ScreenTest for an example of
% how to do detailed checking.

% Open a window and paint the background black
[win, winrect] = Screen('OpenWindow', screenNumber, [0 0 0]*255);
WinTabMex(0, win); %Initialize tablet driver, connect it to 'win'
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %Enable blending to overlay circulor cursor on top of targets
Screen('TextSize',win,30);

MAX_SAMPLES=7e5; %about 20 minutes @ 100Hz = 20*60*100
timevec = zeros(MAX_SAMPLES,1);
thePoints=zeros(MAX_SAMPLES,2);
tabletPoints=zeros(MAX_SAMPLES*2,9); %double the samples since the tablet is sampled @ 200Hz
tabletTime=uint32(zeros(MAX_SAMPLES*2,1)); %collecting from tablet to compare to theoretical times
tablet_queue_length=zeros(MAX_SAMPLES,1);
total_vel=zeros(MAX_SAMPLES,1);
total_displacement=zeros(MAX_SAMPLES,1);
deltax=zeros(MAX_SAMPLES,1);
deltay=zeros(MAX_SAMPLES,1);
time_stamp = NaN(MAX_SAMPLES,8);
timeincirc=zeros(size(psychtest,1),8);
allmaxvel=nan(size(psychtest,1),1);
resx=winrect(3);
resy=winrect(4);
xc= resx-950; %3*resx/4;      % center pixel in x, show on the first monitor % Jordan changed this
yc= resy-600; %resy/2;      % center pixel in y % Jordan changed this
tablet_x_scale = 1/25.4;
tablet_x_offset = -19.2*2540;
tablet_y_scale = -1/25.4;
tablet_y_offset = 12*2540;
%startcol=[160 32 240]; %purple
startcol=[0 255 0]; %we'll make this the same as the wait color (LA)
screenoffset=[resx resy resx resy];
HideCursor; % Hide the mouse cursor.
SetMouse(xc,yc);
screenoffset2=repmat(screenoffset(1:2),4,1); %used to offset the obstacle which has a different size

%%% Features of target file %%%
philist=psychtest(:,1);
x1list=resx-psychtest(:,3);
y1list=resy-psychtest(:,4);

x2list=resx-psychtest(:,9); % 9th and 10th column will be x and y positions of the second target (if we're showing it)
y2list=resy-psychtest(:,10); %add 950 and 800 to the x and y columns for testing

trial_type=psychtest(:,11); %Will range 1-4: (NOTE: ONLY 1 & 3 CURRENTLY USED)
%%%% 1 -> 1 target, no obstacle
%%%% 2 -> 1 target, with obstacle
%%%% 3 -> 2 targets, no obstacle
%%%% 4 -> 2 targets, with obstacle

cued_target=psychtest(:,12); %This will tell us which target to cue on 2 target trials. 0 will signify a leftward target (-30) and 1 will mean rightward (+30).
%Note: just flip a coin as in the flanagan paper?

obstacle_type=psychtest(:,13); %This will tell us what type of barrier to display on a given trial. Take the values of 1, 2 or 3 for each obstacle type. 0 will mean dont show an obstacle.
obs_posdist=psychtest(:,14); %The distribution of which target will be blocked by the obstacle on 2-target trials. 0 will mean block the left target, 1 will mean block the right target
%NOTE: we might want to control this distribution to some degree 
obs_scale = psychtest(:,15); %how much to scale the obstacle length by

vis_feedback=psychtest(:,7);
movement_type=psychtest(:,8);
maxtrialnum=max(movement_type(movement_type<90));
ii=2; %where you want to move to on first trial
hitcirc_count=1;
[~, ~, buttons]=GetMouse;
timeinside=.3; % time in the end target at below vel threshold until mov ends
target_time=0.8; % time of movement (includes timeinside); actual movement time is target_time-timeinside. %%Note: Changed from 0.65 (LA)
wait_time=1; % time in start circle before target appears. Note: Changed from 0.3 (LA)
hold_time=0.75; % time in which experimental conditions of the current trial are displayed before subject initiates movement (LA)
pause_period=wait_time+hold_time; %total time a subject must wait before moving after entering the start target
RT_limit=0.425; %window after completion of ting sound during which subjects must initiate movement (or else a glass sound plays) (LA)
%Generate circle
circlesize=39.37;  %radius of target circle (want 2 cm diameter) (changed from 30-LA)
reward_scale=2; % scaling factor for reward circle around the end point
startcirclesize=circlesize;  %radius of start circle (changed from 15-LA)
cursorsize=circlesize/2;
velthresh=118.11; % velocity theshold (3 cm/sec)
maxtrialtime=5; % 5 seconds;
viewradstart=60; % radius of circle around start that you can see the cursor
dfromstart=200; % how far away you need to be from start to end trial
circletext=Screen('MakeTexture', win, 255*circle(circlesize)); %the target texture
cursortext=Screen('MakeTexture', win, 255*circle(cursorsize));
xstart=x1list(1);       % x and y locations of starting circle
ystart=y1list(1);

cue_target_thresh=78.74; %2 cm in pixels

obs_width=39.37; %1 cm in pixels

xs=resx-xstart; %2nd monitor start positions
ys=resy-ystart;

obs_circletext=Screen('MakeTexture', win, 255/2*circle(obs_width/2)); %circle part of the barrier

% obs_recttext=Screen('MakeTexture', win, 255/2*ones(1,1)); %circle part of the barrier
% rect_loc=[obs_centx-obs_width obs_centy-obs_width obs_centx+obs_width obs_centy+obs_width];

startloc=[xstart-startcirclesize ystart-startcirclesize xstart+startcirclesize ystart+startcirclesize];

insidecircle=0;
started=0;
ismoving=0;
firsttimeout=0;
nextcol1 = [30 144 255]; %light blue
% nextcol2 = [30 144 255];
nextcol2=nextcol1;
tgt_colors=[nextcol1',nextcol2'];
waitcol = [0 255 0];
k = 15;
tab_k = 15;

cue_tgt_flag=0; %the flag that will indicate if the subject has satisfied the criteria to display the cued target.
Score=0;
textxstart_s=resx/2-1850;
textystart_s=5;
obs_check=1; %flag for obstacle feedback
RT_check=1;
%FB_flag=0;
ting_on=1;
xFB_offset=350;
pen_size=5;
MT_flag=0;
%FB_text_shiftx=200;

%%
%%%%%%%%%%%%%%%%%%%%% SETTING UP TRIAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

xloc1_all=0*x1list; % x and y locations of finishing circle
yloc1_all=0*x1list;
xloc2_all=0*x1list;
yloc2_all=0*x1list;


centercircle_x1_all=0*x1list;
centercircle_y1_all=0*x1list;
centercircle_x2_all=0*x1list;
centercircle_y2_all=0*x1list;

p1loc_all=nan(length(x1list),4);
p2loc_all=p1loc_all;
cued_tgt_loc_all=p1loc_all;
ploc_true_all=p1loc_all;

xloc_cue_all=0*x1list;
yloc_cue_all=0*x1list;
centercircle_xcue_all=0*x1list;
centercircle_ycue_all=0*x1list;

xloc_true_all=0*x1list;
yloc_true_all=0*x1list;
centercircle_xtrue_all=0*x1list;
centercircle_ytrue_all=0*x1list;

cued_color_all=zeros(length(x1list),3);
true_color_all=cued_color_all;

obs_dir_all=x1list*0;
current_dir_all=x1list*0;

xloc_obs_all=0*x1list; %the obstacle needs to be positioned indepedently of the target that will be cued
yloc_obs_all=0*x1list;
obs_all=cell(length(x1list),1);
current_angle_all=0*x1list; %we may have cases where we have different angles in 1 block
distp=39.37*20; %these are 20 cm movements

for q=ii:length(x1list)
    xloc1_all(q)=x1list(q);
    yloc1_all(q)= y1list(q);
    xloc2_all(q)=x2list(q);
    yloc2_all(q)=y2list(q);
    
    centercircle_x1_all(q)=xloc1_all(q);
    centercircle_y1_all(q)=yloc1_all(q);
    centercircle_x2_all(q)=xloc2_all(q);
    centercircle_y2_all(q)=yloc2_all(q);
    
    p1loc_all(q,:)=[centercircle_x1_all(q)-circlesize centercircle_y1_all(q)-circlesize centercircle_x1_all(q)+circlesize centercircle_y1_all(q)+circlesize];
    p2loc_all(q,:)=[centercircle_x2_all(q)-circlesize centercircle_y2_all(q)-circlesize centercircle_x2_all(q)+circlesize centercircle_y2_all(q)+circlesize];

    %if its a 2 target trial, determine which target to cue
    if cued_target(q)==0
        cued_tgt_loc_all(q,:)=p1loc_all(q,:);
        cued_color_all(q,:)=nextcol1;
        
        xloc_cue_all(q)=xloc1_all(q);
        yloc_cue_all(q)=yloc1_all(q);
        centercircle_xcue_all(q)=centercircle_x1_all(q);
        centercircle_ycue_all(q)=centercircle_y1_all(q);
        
    elseif cued_target(q)==1
        cued_tgt_loc_all(q,:)=p2loc_all(q,:);
        cued_color_all(q,:)=nextcol2;
        
        xloc_cue_all(q)=xloc2_all(q);
        yloc_cue_all(q)=yloc2_all(q);
        centercircle_xcue_all(q)=centercircle_x2_all(q);
        centercircle_ycue_all(q)=centercircle_y2_all(q);
    end
    
    %determine what the true (absolute) target is depending on the trial type
    if trial_type(q)==1 || trial_type(q)==2
        xloc_true_all(q)=xloc1_all(q);
        yloc_true_all(q)=yloc1_all(q);
        centercircle_xtrue_all(q)=centercircle_x1_all(q);
        centercircle_ytrue_all(q)=centercircle_y1_all(q);
        ploc_true_all(q,:)=[centercircle_xtrue_all(q)-circlesize centercircle_ytrue_all(q)-circlesize centercircle_xtrue_all(q)+circlesize centercircle_ytrue_all(q)+circlesize];
        true_color_all(q,:)=nextcol1;

    elseif trial_type(q)==3 || trial_type(q)==4
        xloc_true_all(q)=xloc_cue_all(q);
        yloc_true_all(q)=yloc_cue_all(q);
        centercircle_xtrue_all(q)=centercircle_xcue_all(q);
        centercircle_ytrue_all(q)=centercircle_ycue_all(q);
        ploc_true_all(q,:)=[centercircle_xtrue_all(q)-circlesize centercircle_ytrue_all(q)-circlesize centercircle_xtrue_all(q)+circlesize centercircle_ytrue_all(q)+circlesize];
        true_color_all(q,:)=cued_color_all(q,:);
    end
    
    %determine if the current trial is a forward or backward one
    if xloc_true_all(q)==xstart & yloc_true_all(q)==ystart
        current_dir_all(q)=0;
    else current_dir_all(q)=1; %1 means its forward
    end
        
     %%%%%%%%%%%% MAKING THE OBSTACLE %%%%%%%%%%%%%%%
     
     %NOTE: should probably move all this stuff to tgt files...
     
     %Need to add offset to obstacles depending on if its a rightward or leftward target
     %We have obstacle A, B and C for each (leftward or rightward) target. B and C have 2 permutationss, where they are flipped against
     %the axis that is the straight line from the start target to the end target
     %Obstacles for straight ahead trials will be handled seprately after the current loop
     
     obs_all{q}.coord=NaN;
     obs_all{q}.edges=NaN;
     
     
     if (current_dir_all(q)==1 && obstacle_type(q)>0), temp_angle=round(abs(asind((resx-xloc_true_all(q)-xs)/distp))); current_angle_all(q)=real(temp_angle); end
     
     if current_angle_all(q) ~=0 %take care of left and right targets first
         if current_angle_all(q)==45
             obs_angle_right=current_angle_all(q)+270;
             obs_angle_left=current_angle_all(q)+0;
         elseif current_angle_all(q)==30
             obs_angle_right=current_angle_all(q)+270;
             obs_angle_left=current_angle_all(q)+30;
         end
         
         if trial_type(q)==3
             %choose which target to position obstacle against
             if obs_posdist(q)==0
                 xloc_obs_all(q)=xloc1_all(q); %left
                 yloc_obs_all(q)=yloc1_all(q);
             elseif obs_posdist(q)==1
                 xloc_obs_all(q)=xloc2_all(q);
                 yloc_obs_all(q)=yloc2_all(q);
             end
         else %position the obs against the target in STT
             xloc_obs_all(q)=xloc_true_all(q);
             yloc_obs_all(q)=yloc_true_all(q);
         end
         
         %Determine the direction of the obstacle (this will help with how we offset the obstacle)
         if (resx-xloc_obs_all(q)-xs)>0
             obs_dir_all(q)=1; %1 means a rightward obs
         elseif (resx-xloc_obs_all(q)-xs)<0
             obs_dir_all(q)=2; %2 means a leftward one
         elseif (resx-xloc_obs_all(q)-xs)==0 && resy-yloc_obs_all(q)~=ys %not the start target
             obs_dir_all(q)=3; %3 means a straight ahead target
         end
         
         if obstacle_type(q)>0 %Note that this should be false if tgt_dir==3
             obs_centx=((resx-xloc_obs_all(q)-xs)/2) + xs; %barrier is centered at the midpoint between a target and the start position
             obs_centy=((resy-yloc_obs_all(q)-ys)/2) + ys;
             %xxx(q)=obs_centx; %checking if the obstacle is centered correctly
             %yyy(q)=obs_centy;
             obs_length = 39.37 * obs_scale(q);
             
             if obstacle_type(q)==1 %Obstacle A
                 
                 if current_angle_all(q)==45
                     obs_length1=39.37*9; %modulate length for obs 1 based on angle of movement
                 else
                     obs_length1=39.37*14;
                 end
                 
                 if obs_dir_all(q)==1 %right
                     [P1temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_right,obs_width,obs_length1);
                     edge_point=(P1temp(1,:) + P1temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,2*39.37,obs_length1,obs_dir_all(q),1);
                     [P1,~,~]=make_obs_loc(ncx,ncy,obs_angle_right,obs_width,obs_length1);
                     obs_all{q}.coord=P1;
                     obs_all{q}.edges={[P1(3,:);P1(2,:)],[P1(4,:);P1(1,:)],[P1(3,:);P1(4,:)],[P1(2,:);P1(1,:)]};
                     
                 elseif obs_dir_all(q)==2 %left
                     [P1temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_left,obs_width,obs_length1);
                     edge_point=(P1temp(1,:) + P1temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,2*39.37,obs_length1,obs_dir_all(q),1);
                     [P1,~,~]=make_obs_loc(ncx,ncy,obs_angle_left,obs_width,obs_length1);
                     obs_all{q}.coord=P1;
                     obs_all{q}.edges={[P1(3,:);P1(2,:)],[P1(4,:);P1(1,:)],[P1(3,:);P1(4,:)],[P1(2,:);P1(1,:)]};
                     
                 end
                 
                 %Circle stuff
                 %C1_loc=[C1(1)-obs_width C1(2)-obs_width C1(1)+obs_width C1(2)+obs_width];
                 %C1_flip_loc=[C1_flip(1)-obs_width C1_flip(2)-obs_width C1_flip(1)+obs_width C1_flip(2)+obs_width];
                 
                 
                 %             obs_circ_loc=C1_loc;
                 %             obs_circ_flip_loc=C1_flip_loc;
                 %             obs_circx=C1(1);
                 %             obs_circy=C1(2);
                 
             elseif obstacle_type(q)==2 || obstacle_type(q)==3 %Obstacle B
                 %if q==4,keyboard;end
                 obs_length2=obs_length; %2 cm long, so shift by half the length, so enter 0 for 3rd argument
                 
                 if obs_dir_all(q)==1 && obstacle_type(q)==2 %right
                     [P2temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_right,obs_width,obs_length2);
                     edge_point=(P2temp(1,:) + P2temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,0,obs_length2,obs_dir_all(q),2);
                     [P2,~,~]=make_obs_loc(ncx,ncy,obs_angle_right,obs_width,obs_length2);
                     obs_all{q}.coord=P2;
                     obs_all{q}.edges={[P2(3,:);P2(2,:)],[P2(4,:);P2(1,:)],[P2(3,:);P2(4,:)],[P2(2,:);P2(1,:)]};
                     
                 elseif obs_dir_all(q)==2 && obstacle_type(q)==2
                     [P2temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_left,obs_width,obs_length2);
                     edge_point=(P2temp(1,:) + P2temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,0,obs_length2,obs_dir_all(q),2);
                     [P2,~,~]=make_obs_loc(ncx,ncy,obs_angle_left,obs_width,obs_length2);
                     obs_all{q}.coord=P2;
                     obs_all{q}.edges={[P2(3,:);P2(2,:)],[P2(4,:);P2(1,:)],[P2(3,:);P2(4,:)],[P2(2,:);P2(1,:)]};
                     
                 elseif obs_dir_all(q)==1 && obstacle_type(q)==3 %symmetric obstacle type
                     [P2temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_right,obs_width,obs_length2);
                     edge_point=(P2temp(1,:) + P2temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,0,obs_length2,obs_dir_all(q),1);
                     [P2,~,~]=make_obs_loc(ncx,ncy,obs_angle_right,obs_width,obs_length2);
                     obs_all{q}.coord=P2;
                     obs_all{q}.edges={[P2(3,:);P2(2,:)],[P2(4,:);P2(1,:)],[P2(3,:);P2(4,:)],[P2(2,:);P2(1,:)]};
                     
                 elseif obs_dir_all(q)==2 && obstacle_type(q)==3
                     [P2temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_left,obs_width,obs_length2);
                     edge_point=(P2temp(1,:) + P2temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,0,obs_length2,obs_dir_all(q),1);
                     [P2,~,~]=make_obs_loc(ncx,ncy,obs_angle_left,obs_width,obs_length2);
                     obs_all{q}.coord=P2;
                     obs_all{q}.edges={[P2(3,:);P2(2,:)],[P2(4,:);P2(1,:)],[P2(3,:);P2(4,:)],[P2(2,:);P2(1,:)]};
                 end
                 
                 %C2_loc=[C2(1)-obs_width C2(2)-obs_width C2(1)+obs_width C2(2)+obs_width];
                 %C2_flip_loc=[C2_flip(1)-obs_width C2_flip(2)-obs_width C2_flip(1)+obs_width C2_flip(2)+obs_width];
                 
                 
                 %             obs_circ_loc=C2_loc;
                 %             obs_circ_flip_loc=C2_flip_loc;
                 %             obs_circx=C2(1);
                 %             obs_circy=C2(2);
                 
             elseif obstacle_type(q)==4 || obstacle_type(q)==5 %Obstacle C
                 
                 obs_length3=obs_length; %a square
                 
                 if obs_dir_all(q)==1 && obstacle_type(q)==4 %right
                     %shift by -1/2L since we want 1 cm between the obstacle and the axis
                     [P3temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_right,obs_width,obs_length3);
                     edge_point=(P3temp(1,:) + P3temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,-(obs_length3),obs_length3,obs_dir_all(q),2);
                     [P3,~,~]=make_obs_loc(ncx,ncy,obs_angle_right,obs_width,obs_length3);
                     obs_all{q}.coord=P3;
                     obs_all{q}.edges={[P3(3,:);P3(2,:)],[P3(4,:);P3(1,:)],[P3(3,:);P3(4,:)],[P3(2,:);P3(1,:)]};
                     
                 elseif obs_dir_all(q)==2 && obstacle_type(q)==4 %left
                     [P3temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_left,obs_width,obs_length3);
                     edge_point=(P3temp(1,:) + P3temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,-(obs_length3),obs_length3,obs_dir_all(q),2);
                     [P3,~,~]=make_obs_loc(ncx,ncy,obs_angle_left,obs_width,obs_length3);
                     obs_all{q}.coord=P3;
                     obs_all{q}.edges={[P3(3,:);P3(2,:)],[P3(4,:);P3(1,:)],[P3(3,:);P3(4,:)],[P3(2,:);P3(1,:)]};
                     
                 elseif obs_dir_all(q)==1 && obstacle_type(q)==5
                     [P3temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_right,obs_width,obs_length3);
                     edge_point=(P3temp(1,:) + P3temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,-(obs_length3),obs_length3,obs_dir_all(q),1);
                     [P3,~,~]=make_obs_loc(ncx,ncy,obs_angle_right,obs_width,obs_length3);
                     obs_all{q}.coord=P3;
                     obs_all{q}.edges={[P3(3,:);P3(2,:)],[P3(4,:);P3(1,:)],[P3(3,:);P3(4,:)],[P3(2,:);P3(1,:)]};
                     
                 elseif obs_dir_all(q)==2 && obstacle_type(q)==5
                     [P3temp,~,~]=make_obs_loc(obs_centx,obs_centy,obs_angle_left,obs_width,obs_length3);
                     edge_point=(P3temp(1,:) + P3temp(2,:))/2;
                     [ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,-(obs_length3),obs_length3,obs_dir_all(q),1);
                     [P3,~,~]=make_obs_loc(ncx,ncy,obs_angle_left,obs_width,obs_length3);
                     obs_all{q}.coord=P3;
                     obs_all{q}.edges={[P3(3,:);P3(2,:)],[P3(4,:);P3(1,:)],[P3(3,:);P3(4,:)],[P3(2,:);P3(1,:)]};
                 end
                 
                 %C3_loc=[C3(1)-obs_width C3(2)-obs_width C3(1)+obs_width C3(2)+obs_width];
                 %C3_flip_loc=[C3_flip(1)-obs_width C3_flip(2)-obs_width C3_flip(1)+obs_width C3_flip(2)+obs_width];
                 
                 
                 %             obs_circ_loc=C3_loc;
                 %             obs_circ_flip_loc=C3_flip_loc;
                 %             obs_circx=C3(1);
                 %             obs_circy=C3(2);
             end
         end
     end
end

%find each obstacle type to reuse for the straight ahead case
obs_30_left_idx=find(obstacle_type>0 & current_angle_all==30 & obs_dir_all==2 & trial_type==1,1,'first');
obs_30_right_idx=find(obstacle_type>0 & current_angle_all==30 & obs_dir_all==1 & trial_type==1,1,'first');
obs_45_left_idx=find(obstacle_type>0 & current_angle_all==45 & obs_dir_all==2 & trial_type==1,1,'first');
obs_45_right_idx=find(obstacle_type>0 & current_angle_all==45 & obs_dir_all==1 & trial_type==1,1,'first');

if isempty(obs_30_left_idx) & isempty(obs_30_right_idx) & isempty(obs_45_left_idx) & isempty(obs_45_right_idx),
    obs_block=0; %then we're not in a block w/ obs trials
else
    obs_block=1;
end

if obs_block,
    if isempty(obs_45_left_idx),
        is45=0;
        obs_30_left_constr=obs_all{obs_30_left_idx};
        obs_30_right_constr=obs_all{obs_30_right_idx};
    else
        is45=1;
        obs_45_left_constr=obs_all{obs_45_left_idx};
        obs_45_right_constr=obs_all{obs_45_right_idx};
    end
    
    for qw=ii:length(x1list)
        if current_angle_all(qw)==0 && obstacle_type(qw)>0
            %we will determine where to place the obstacle for these trials based on obs_posdist (14th column)
            if is45==0 %then its a 30 degree block
                if obs_posdist(qw)==0
                    obs_all{qw}=obs_30_left_constr;
                elseif obs_posdist(qw)==1
                    obs_all{qw}=obs_30_right_constr;
                end
            elseif is45==1
                if obs_posdist(qw)==0
                    obs_all{qw}=obs_45_left_constr;
                elseif obs_posdist(qw)==1
                    obs_all{qw}=obs_45_right_constr;
                end
            end
        end
    end
end
obs_coord=obs_all{ii}.coord;
xloc_true=xloc_true_all(ii);
yloc_true=yloc_true_all(ii);
p1loc=p1loc_all(ii,:);
p2loc=p2loc_all(ii,:);
ploc_true=ploc_true_all(ii,:);
% cued_tgt_loc=cued_tgt_loc_all(ii,:);
% cued_color=cued_color_all(ii,:);
centercircle_xtrue=centercircle_xtrue_all(ii);
centercircle_ytrue=centercircle_ytrue_all(ii);
true_color=true_color_all(ii,:);
all_tgt_loc=[p1loc',p2loc'];
current_dir=current_dir_all(ii);

bad_label=zeros(length(x1list),3); %flag if the subject jumped the gun, initiated the movement too slowly, or hit the obstacle
MT_label=zeros(length(x1list),1); %flag if the subject got the movement with the correct movement time
tgt_leave_time=nan(length(x1list),1); %record the times where the subject first leaves the target
tgt_leave_flag=1;
ting_time_onset=nan(length(x1list),1); %we record the time at which the ting goes off to calculate reaction time
movement_onset_time=nan(length(x1list),1); %Need to record the movement onset times to calculate reaction time

%% Main loop 
WinTabMex(2); %Empties the packet queue in preparation for collecting actual data
% WaitSecs(0.25); % just in case wait
startTime = GetSecs;
t_a=startTime;
movement_start_time=GetSecs;
tic;
while ~any(buttons(2:end))&& ii<length(philist)
    %ddd=[xxx(ii)-8 yyy(ii)-8 xxx(ii)+8 yyy(ii)+8];
    k=k+1;
    time_stamp(k,1) = toc;
    Screen('Flip',win,0,0,2);
    
%     if ii==8 && started==1
%         keyboard;
%     end

%     if started==1, keyboard; end;
    
    %%% Collect tablet samples %%%
%    [mX, mY, buttons] = GetMouse; old code 
    [~, ~, buttons] = GetMouse;
    time_stamp(k,2) = toc;  
    pkt = WinTabMex(5);
    tablet_queue_length(k) = 0;
    
    while ~isempty(pkt)
        tabletPoints(tab_k,1:8) = pkt(1:8)';
        tabletPoints(tab_k,9) = toc;
        tabletTime(tab_k) = pkt(6)-tabletTime(16); %in ms
        tab_k = tab_k+1;
        tablet_queue_length(k) = tablet_queue_length(k)+1;
        pkt = WinTabMex(5);
    end
    mX = (tabletPoints(tab_k-1,1)-tablet_x_offset)*tablet_x_scale;
    mY = (tabletPoints(tab_k-1,2)-tablet_y_offset)*tablet_y_scale;
        
    time_stamp(k,3) = toc;
    thePoints(k,:) = [mX mY]; % record full precision points
    t=GetSecs;
    timevec(k)=t-startTime;
    % location to display cursor
    [dx dy]=rotatexy(round(mX)-xc,-(round(mY)-yc),philist(ii));         % rotate xy centers around middle of screen coords at (0,0) and round to screen precision
    dx = dx+xc;     % shift screen back
    dy = -dy+yc;     % shift screen back
    % Draw the sprite at the new location.
    cursor=[min(max(dx-cursorsize,resx/2),resx-10) min(max(dy-cursorsize,0),resy-10) max(min(dx+cursorsize,resx),resx/2+10) max(min(dy+cursorsize,resy),10)];
    time_stamp(k,4) = toc;
         
    %%%%%%%%%%%%% MOVEMENT STUFF %%%%%%%%%%%%%%%%%%
    
    if started==0    % inside start circle wait time (not yet started)
        
        if movement_type(ii) ~= 99
            Screen('DrawTexture', win, circletext,[],startloc,[],[],[],waitcol); %Drawing the target at the start location for the subject screen
            Screen('DrawTexture', win, circletext,[],screenoffset-startloc(1,[3,4,1,2]),[],[],[],waitcol); %Drawing the target at the start location for the experimenter screen
        end
        if sqrt((dx-xstart)^2+(dy-ystart)^2)>circlesize && current_dir==1     % cursor left old target circle, removed reward scale (LA - 8/4/16). NOTE: Should this be startcirclesize?
            if firsttimeout==1 && (movement_type(ii)~=99)
                firsttimeout=0;
                Snd('Play',glass,fs1);
                %FB_flag=1;
                ting_on=0;
                bad_label(ii,1)=1;
            end
            movement_start_time=GetSecs;
            t_a=GetSecs;
        elseif sqrt((dx-xstart)^2+(dy-ystart)^2)>startcirclesize
            movement_start_time=GetSecs;
            t_a=GetSecs;
        end
        
        if (GetSecs >= t_a + wait_time)
            
            Screen('DrawTexture', win, circletext,[],startloc,[],[],[],startcol); 
            Screen('DrawTexture', win, circletext,[],screenoffset-startloc(1,[3,4,1,2]),[],[],[],startcol);
            
            if obstacle_type(ii)>0
                Screen('FillPoly', win,grey, obs_coord,1); % convex polygon
                Screen('FillPoly', win,grey, screenoffset2-[obs_coord(3,:);obs_coord(4,:);obs_coord(1,:);obs_coord(2,:)],1);
                %Screen('DrawTexture', win, circletext,[],ddd,[],[],[],[255,0,0]);
                %Screen('DrawTexture', win, circletext,[],screenoffset-ddd(1,[3,4,1,2]),[],[],[],[255,0,0]);
                %Screen('DrawTexture', win,obs_circletext,[],obs_circ_loc,[],[],[],grey);
                %Screen('DrawTexture', win,obs_circletext,[],screenoffset-obs_circ_loc(1,[3,4,1,2]),[],[],[],grey);
            end
            
            if trial_type(ii)==1 || trial_type(ii)==2 %target should be cued (filled in) immediately if its a single target trial
                Screen('DrawTexture', win,circletext,[],ploc_true,[],[],[],true_color);
                Screen('DrawTexture', win, circletext,[],screenoffset-ploc_true(1,[3,4,1,2]),[],[],[],true_color);
                
            elseif trial_type(ii) ==3 || trial_type(ii)==4
                Screen('FrameOval',win,tgt_colors,all_tgt_loc,pen_size); %both targets
                Screen('FrameOval',win,nextcol1,screenoffset-p1loc(1,[3,4,1,2]),pen_size);
                Screen('FrameOval',win,nextcol2,screenoffset-p2loc(1,[3,4,1,2]),pen_size);
            end   
        end
        %flag if the participant has started depending on if its a forward/backward trial
        if current_dir==0
            ting_on=0; %dont play feedback if its a backward trial
            %Snd('Play',ting,fs);
            firsttimeout=1;
            started=1;
            t_a=GetSecs;
        elseif (current_dir==1) && ((GetSecs >= t_a + pause_period) || movement_type(ii)==99)  % in circle for long enough
            Snd('Play',ting,fs);
            ting_time_onset(ii)=GetSecs;
            firsttimeout=1;
            started=1;
            t_a=GetSecs;
        end
        
    else    %  cumulative wait time is over
        %keyboard;
        %check if you started too slow
        %IMPORTANT: The RT limit has been referenced to movement onset (i.e. 2 cm into the movement)-this may be too harsh - 8/8/16
        if (GetSecs>=(t_a+RT_limit+ting_time)) && (sqrt((dx-xstart)^2+(dy-ystart)^2)<(2*circlesize)) && (RT_check==1) && current_dir==1 && cue_tgt_flag==0
            Snd('Play',car_alarm,fs5);
            RT_check=0;
            %FB_flag=2;
            bad_label(ii,2)=1;
            ting_on=0;
            firsttimeout=0;
            started=0;
        end
        
        if obstacle_type(ii)>0
            Screen('FillPoly', win,grey, obs_coord,1); % convex polygon
            Screen('FillPoly', win,grey, screenoffset2-[obs_coord(3,:);obs_coord(4,:);obs_coord(1,:);obs_coord(2,:)],1);
            %Screen('DrawTexture', win, circletext,[],ddd,[],[],[],[255,0,0]);
            %Screen('DrawTexture', win, circletext,[],screenoffset-ddd(1,[3,4,1,2]),[],[],[],[255,0,0]);
            %Screen('DrawTexture', win,obs_circletext,[],obs_circ_loc,[],[],[],grey);
            %Screen('DrawTexture', win,obs_circletext,[],screenoffset-obs_circ_loc(1,[3,4,1,2]),[],[],[],grey);
        end
                
        if trial_type(ii)==1 || trial_type(ii)==2
            Screen('DrawTexture', win,circletext,[],ploc_true,[],[],[],true_color);
            Screen('DrawTexture', win, circletext,[],screenoffset-ploc_true(1,[3,4,1,2]),[],[],[],true_color);
            
            %Screen('DrawLine',win,nextcol1,xstart,ystart,xloc_true,yloc_true); %%%CHECKING OBSTACLE POSITIONING%%%
            
        elseif (trial_type(ii)==3 || trial_type(ii)==4) && cue_tgt_flag==0
            Screen('FrameOval',win,tgt_colors,all_tgt_loc,pen_size); %both targets
            Screen('FrameOval',win,nextcol1,screenoffset-p1loc(1,[3,4,1,2]),pen_size);
            Screen('FrameOval',win,nextcol2,screenoffset-p2loc(1,[3,4,1,2]),pen_size);
            
            %Screen('DrawLine',win,nextcol1,xstart,ystart,xloc1_all(ii),yloc1_all(ii)); %%%CHECKING OBSTACLE POSITIONING%%%
            %Screen('DrawLine',win,nextcol1,xstart,ystart,xloc2_all(ii),yloc2_all(ii)); %%%CHECKING OBSTACLE POSITIONING%%%
        end
        
        if cue_tgt_flag %cue needs to happen on all trials
            Screen('DrawTexture', win,circletext,[],ploc_true,[],[],[],true_color);
            Screen('DrawTexture', win, circletext,[],screenoffset-ploc_true(1,[3,4,1,2]),[],[],[],true_color);
            %                 if ii==8,keyboard; end
        end
        
        if movement_type(ii) ~= 99 % make purple circle at start location
            Screen('DrawTexture', win, circletext,[],startloc,[],[],[],startcol);
            Screen('DrawTexture', win, circletext,[],screenoffset-startloc(1,[3,4,1,2]),[],[],[],startcol);
        end
        %Check if it is inside target circle
        if sqrt((dx-centercircle_xtrue)^2+(dy-centercircle_ytrue)^2)<circlesize && insidecircle==0 && ismoving==1    % first time in circle
            t1=GetSecs;
            insidecircle=1; %circlesize was multiplied reward_scale for insidecircle (LA changed this)
        elseif sqrt((dx-centercircle_xtrue)^2+(dy-centercircle_ytrue)^2)<circlesize && ismoving==1 
            insidecircle=1; %NOTE: may want to introduce a reward scale for forward vs backward instead of getting rid of it altogether (LA)
        else
            insidecircle=0;
        end
        tmp_GS = GetSecs; %When used, a clock counting the trial duration      
        if ((insidecircle==1) && (tmp_GS-t1 >= timeinside) && (ismoving == 1)) ... %Note: Maybe make it so that you're inside the circle long engouh at below the vel threshold? (LA)
                || ((~vis_feedback(ii) && movement_type(ii) < 99) ...
                && (tmp_GS-movement_start_time > maxtrialtime ...
                || (max(total_vel((k-12):k)) < velthresh && sqrt((dx-xc)^2 + (dy-yc)^2) > dfromstart)))
            Snd('Wait');    % (inside circle for long enough) or (no feedback and (took too long or (traveled long enough and stopped for long enough)))
            Snd('Quiet');   % move on to next movement
            Snd('Close');
            
            if (insidecircle==1) && (hitcirc_count>20)
                
                tmp_c = 1;
                tmp_a = tmp_GS-movement_start_time;
                tmp_b = mean(timeincirc(hitcirc_count-(2:2:20),4));
                
            elseif (insidecircle==1) && (hitcirc_count>2)
                tmp_c = 2;
                tmp_a = tmp_GS-movement_start_time;
                tmp_b = mean(timeincirc(hitcirc_count-2:-2:1,4));
                
            elseif  (insidecircle==1) && (hitcirc_count<=2)
                tmp_c = 3;
                tmp_a = NaN;
                tmp_b = NaN;
            end
            if (tmp_GS-movement_start_time<=(target_time)) && ting_on
                res=0;
                nextcol1 = [0 0 255];
                Snd('Play',mario_coin,fs2);
                Score=Score+1;
                MT_label(ii)=1;
            else
                res=1;
                %nextcol1 = [255 0 0];
            end
            if current_dir==1 & (tmp_GS-movement_start_time>(target_time))
                MT_flag=1;
                Snd('Play',MT_bad_buzzer,fs6);
                MT_FBx=xloc_true;
                MT_FBy=yloc_true;
            end
           
            timeincirc(hitcirc_count,:)=[movement_start_time-startTime movement_start_time-startTime tmp_GS-startTime tmp_GS-movement_start_time res tmp_c tmp_a tmp_b];
            hitcirc_count=hitcirc_count+1;
            ii=ii+1; %current target # / movement #
            xstart=xloc_true;                     % coords for new oval(s)
            ystart=yloc_true;
            
            startloc=[xstart-startcirclesize ystart-startcirclesize xstart+startcirclesize ystart+startcirclesize];
            %keyboard;
            obs_coord=obs_all{ii}.coord;
            xloc_true=xloc_true_all(ii);
            yloc_true=yloc_true_all(ii);
            p1loc=p1loc_all(ii,:);
            p2loc=p2loc_all(ii,:);
            ploc_true=ploc_true_all(ii,:);
            true_color=true_color_all(ii,:);
            %cued_tgt_loc=cued_tgt_loc_all(ii,:);
            %cued_color=cued_color_all(ii,:);
            centercircle_xtrue=centercircle_xtrue_all(ii);
            centercircle_ytrue=centercircle_ytrue_all(ii);
            all_tgt_loc=[p1loc',p2loc'];
            current_dir=current_dir_all(ii);
            
            t_a=tmp_GS;
            t1=tmp_GS;
            insidecircle=0;
            started=0;
            ismoving=0;
            cue_tgt_flag=0;
            obs_check=1;
            RT_check=1;
            %FB_flag=0;
            if rem(ii,2)==0, MT_flag=0; end
            ting_on=1; %a flag that prevents feedback from playing if the obstacle is hit or if its a return trial
            tgt_leave_flag=1;
        end
    end
    time_stamp(k,5) = toc;
    % if visual feedback present, draw cursor, otherwise don't
    % if cursor is around viewradstart radius circle of the middle in type 99 movements
    if ~(movement_type(ii)==99 && sqrt((dx-xc)^2 + (dy-yc)^2) > viewradstart) && (~ismoving || vis_feedback(ii) || ((movement_type(ii)==99 || ~started) && sqrt((dx-xc)^2 + (dy-yc)^2) < viewradstart))
        %Screen('DrawTexture', win, cursortext,[],cursor,[],[],[],[255 255 255]);
        %Screen('DrawTexture', win, cursortext,[],screenoffset-cursor(1,[3,4,1,2]),[],[],[],[255 255 255]);
        Screen('DrawDots',win,[resx-dx;resy-dy],cursorsize,[255,255,255],[],1);
        Screen('DrawDots',win,[dx;dy],cursorsize,[255,255,255],[],1);
    else
        Screen('DrawTexture', win, cursortext,[],screenoffset-cursor(1,[3,4,1,2]),[],[],[],[255 0 255]);
    end
    
%     deltaxlast=diff(thePoints(k-1:k,1));
%     deltaylast=diff(thePoints(k-1:k,2));
    
    %NOTE:tablet data used here instead of screen (LA - 8/4/16)
    deltaxlast=diff(double(tabletPoints(tab_k-1:tab_k,1)))*tablet_x_scale; 
    deltaylast=diff(double(tabletPoints(tab_k-1:tab_k,2)))*tablet_y_scale;
    total_displacementlast=sqrt(deltaxlast.^2+deltaylast.^2);
    
    total_vellast=total_displacementlast/diff(timevec(k-1:k))';
    deltax(k)=deltaxlast;
    deltay(k)=deltaylast;
    total_displacement(k)=total_displacementlast;
    total_vel(k)=total_vellast;
   
    % if no longer waiting, movement starts when moving past 2 cm (removed velocity threshold criteria 8/5/16 -LA)
    if started==1 && insidecircle==0 && ismoving == 0 && (sqrt((dx-xstart)^2+(dy-ystart)^2) >= cue_target_thresh)
        movement_start_time=GetSecs;
        ismoving = 1;
        movement_onset_time(ii)=movement_start_time;
    elseif ismoving == 0
        movement_start_time=GetSecs;
    end
    
        %%%%%%%%%%%%% CUE CRITERIA %%%%%%%%%%%%
        %NOTE: Changed it so that we no longer check this criteria once the flag goes off
      if started==1 && sqrt((dx-xstart)^2+(dy-ystart)^2) >= cue_target_thresh && cue_tgt_flag==0 
          cue_tgt_flag=1;  
      end
      
      %%%%%%% TARGET LEAVING TIMES %%%%%%%%%%%%
      if started==1 && insidecircle==0 && tgt_leave_flag && (sqrt((dx-xstart)^2+(dy-ystart)^2) >= startcirclesize)
          tgt_leave_time(ii)=GetSecs;
          tgt_leave_flag=0;
      end
      
      %%%%%%%%%% VELOCITY THRESHOLD TIMES %%%%%%%%%%%%%%
      %%% record the time at which the velocity exceeds a threshold and stays above that threshold for a certain amount of time
      %NOTE: THIS NEEDS DEBUGGING
%       if started==1 && insidecircle==0 && (total_vellast > velthresh) && vel_check==1
%           veltime_flag=GetSecs;
%           tempvel=total_vellast;
%           vel_check=0;
%           veltime_check=1;
%       end
%       while GetSecs-veltime_flag<0.2 & veltime_check==1
%           if total_vellast>=tempvel
%               veltime(ii)=veltime_flag;
%               veltime_check=0;
%           else
%               vel_check=1;  
%           end
%       end
    
    %%%%%%%%%% CHECK IF INSIDE BARRIER %%%%%%%%%%
    if obstacle_type(ii)>0 & started==1, 
        obstacle_flag=detect_collision_V8(obs_all{ii}.coord,[resx-mX,resy-mY],obs_all{ii}.edges); 
    end 

    if (started==1 && obstacle_type(ii)>0 && obstacle_flag && obs_check==1), %need to add the circle part
         Snd('Play',obstacle_error,fs3);
         obs_check=0;
         Score=Score-10;
         %FB_flag=3;
         bad_label(ii,3)=1;
         ting_on=0;
%          keyboard;
    end
    
    %SHOWING SCORE
    Screen('DrawText', win, ['YOUR SCORE: ' num2str(Score)], textxstart_s, textystart_s,[255 20 147]);
    dummy_=DrawFormattedText(win, ['YOUR SCORE: ' num2str(Score)], resx/2+1600, resy/2+400,[255 20 147],[],1,1,[],[]);
    
    % OBSTACLE TEXT FEEDBACK
    if bad_label(ii,3)==1
        Screen('DrawText', win, ['-10!'], textxstart_s+xFB_offset, textystart_s,[255 255 255]);
        dummy2_=DrawFormattedText(win, ['-10!'], resx/2+1600-2*xFB_offset, resy/2+400,[255 255 255],[],1,1,[],[]);
    end
    

    %% MOVEMENT VISUAL FEEDBACK (text) 
    %NOTE: One condition (e.g., too early) doesnt overide the other -> fix this?
%     if (FB_flag>0) && (current_dir==1) && insidecircle==0
%         if FB_flag==1
%             Screen('DrawText', win, ['Too early!'], xs+FB_text_shiftx, ys,[255 255 255]);
%             dummy2_=DrawFormattedText(win, ['Too early!'], resx-xs+FB_text_shiftx, resy-ys,[255 255 255],[],1,1,[],[]);
%         elseif FB_flag==2
%             Screen('DrawText', win, ['Too late!'], xs-FB_text_shiftx-100, ys,[255 255 255]);
%             dummy2_=DrawFormattedText(win, ['Too late!'], resx-xs-FB_text_shiftx-100, resy-ys,[255 255 255],[],1,1,[],[]);
%         end
%     end
    
    if MT_flag==1 %making this separate so it doesnt get overided
        Screen('DrawText', win, ['Too slow!'], resx-MT_FBx-60, resy-MT_FBy-88,[255 255 255]);
        dummy2_=DrawFormattedText(win, ['Too Slow!'], MT_FBx-60, ...
            MT_FBy+78,[255 255 255],[],1,1,[],[]); %exp screen
    end

%%         
    % Draw additional info in the experiment window
    time_stamp(k,6) = toc;
    if (hitcirc_count>1) %%&& (mod(hitcirc_count,2)==0)
        textxstart=resx/2-800;
        textystart=10;
        ydist=20;
        if isempty(find(timevec(1:k)>timeincirc(hitcirc_count-1,3),1))
            t_end=k;
        else
            t_end=find(timevec(1:k)>timeincirc(hitcirc_count-1,3),1)-1;
        end
        t_start=find(timevec(1:k)>timeincirc(hitcirc_count-1,1),1)-1;
        lastmaxvel=max(smooth(total_vel(t_start:t_end)))*.0254/100;
        allmaxvel(hitcirc_count)=lastmaxvel;
        lastpath=repmat([resx;resy],size(t_start:t_end))-thePoints(t_start:t_end,:)';
        lastpath2=lastpath(:,sort([1 2:length(lastpath) 2:length(lastpath)]));
        lastvel=smooth(total_vel(t_start:t_end))';
        velxvals=linspace(0,500,length(lastvel));
        velyvals=-lastvel/max(lastvel)*100;
        lastvel2=[velxvals;velyvals];
        lastvel3=lastvel2(:,sort([1 2:length(lastvel2) 2:length(lastvel2)]));
        Screen('DrawLines', win, lastpath2,1,[255 255 0]);
        Screen('DrawLines', win, lastvel3,1,[255 255 0],[textxstart textystart+10*ydist]);
        Screen('DrawText', win, ['Movetime: ' num2str(timeincirc(hitcirc_count-1,3)-timeincirc(hitcirc_count-1,2))], textxstart, textystart,[255 0 0]);
        Screen('DrawText', win, ['Max Vel: ' num2str(lastmaxvel)], textxstart, textystart+ydist,[255 0 0]);
        Screen('DrawText', win, ['Ave Movetime (block): ' num2str(nanmean(timeincirc(1:2:hitcirc_count-1,3)-timeincirc(1:2:hitcirc_count-1,2)))], textxstart, textystart+2*ydist,[255 0 0]);
        Screen('DrawText', win, ['Ave Max Vel (block): ' num2str(nanmean(allmaxvel(3:end)))], textxstart, textystart+3*ydist,[255 0 0]);
                
    end
    time_stamp(k,7) = toc;
end

%% Clean-up + save
ShowCursor;
Screen('CloseAll'); % close screen
WinTabMex(3); % Stop/Pause data acquisition.
WinTabMex(1); % Shutdown driver.

timevec=timevec(1:k); % difference between timevec and timevec2 is time it takes to run while loop
thePoints=thePoints(1:k,:);
tablet_queue_length=tablet_queue_length(1:k);
tabletPoints=tabletPoints(1:tab_k,:); %#ok<NASGU>
tabletTime=tabletTime(1:tab_k,:); %added 8/8/16 (LA)
time_stamp = time_stamp(1:k,:); %#ok<NASGU>
total_vel=total_vel(1:k); %#ok<NASGU>
total_displacement=total_displacement(1:k); %#ok<NASGU>
deltax=deltax(1:k); %#ok<NASGU>
deltay=deltay(1:k); %#ok<NASGU>
timeincirc=timeincirc(1:hitcirc_count-1,:);
thePoints(:,1)=thePoints(:,1)-resx/2; %#ok<NASGU>        % adjust for other monitor points
goodlist=find(timeincirc(:,5)==0);
typelist=cell(maxtrialnum,1);

for counter=1:maxtrialnum
    typelist{counter}=find(movement_type==counter)-1;
end
name_prefix_all = [name_prefix,'_set_',tgt_set,'_',date];
disp('Saving...')
if ~exist([name_prefix_all,'.mat'],'file'), datafile_name = [name_prefix_all,'.mat'];
elseif ~exist([name_prefix_all,'_a.mat'],'file'), datafile_name = [name_prefix_all,'_a.mat'];
elseif ~exist([name_prefix_all,'_b.mat'],'file'), datafile_name = [name_prefix_all,'_b.mat'];
else
    char1='c';
    while exist([name_prefix_all,'_',char1,'.mat'],'file'), char1=char(char1+1); end
    datafile_name = [name_prefix_all,'_',char1,'.mat'];
end
save(datafile_name); disp(['Saved ', datafile_name]);


disp(['You got a total of ' num2str(length(goodlist)) ' trials correct.'])
curlen=1;
maxlen=0;
for ii=1:length(goodlist)
    if curlen>maxlen
        maxlen=curlen;
    end
    
    if ii>1
        if goodlist(ii)==goodlist(ii-1)+1
            curlen=curlen+1;
        else
            curlen=1;
        end
    end
end
disp(['Your longest streak was ' num2str(maxlen) ' trials in a row.'])
disp(['You took ' num2str(timevec(end)) ' seconds to complete this block.'])
return


%% Extra functions
function [rx, ry] = rotatexy(x,y,phi)
% phi is in degrees
phi=phi*pi/180;
[theta r]=cart2pol(x,y);
[rx ry]=pol2cart(theta+phi,r);
return

function [rect_points,circle_center,circle_center_flip]=make_obs_loc(cx,cy,theta,width,height)
%function to define the coordinates of a rotated rectangle, drawn to the
%screen as a polygon, NOT a texture
%this function will also return a 2nd point known to be at the center of
%one of the 2 short edges
theta=deg2rad(theta);
UL = [cx + (width/2) * cos(theta) - (height/2) * sin(theta) , cy + (height/2) * cos(theta) + (width/2) * sin(theta)];
UR = [cx - (width/2) * cos(theta) - (height/2) * sin(theta) , cy + (height/2) * cos(theta) - (width/2) * sin(theta)];
LL = [cx + (width/2) * cos(theta) + (height/2) * sin(theta) , cy - (height/2) * cos(theta) + (width/2) * sin(theta)];
LR = [cx - (width/2) * cos(theta) + (height/2) * sin(theta) , cy - (height/2) * cos(theta) - (width/2) * sin(theta)];

rect_points=[UL;UR;LR;LL]; %order is important here
%rect_points=[LL;LR;UR;UL];

circle_center=((UR-UL)/2)+UL; %check this
circle_center_flip=((LR-LL)/2) +LL;

return

function [new_cx,new_cy]=calculate_obstacle_offset(cx,cy,Q2temp,P,obs_length,dir,T)
%this function offsets the center of the obstacle so that it is protruding P pixels past the axis that is
%the straight line between the start and end target

%This assumes that cx and cy are passed as the center coordinates between the start and end target (or wherever the center point of the rectangle is)

%We'll do this by obtaining a vector v the starts at the center given, and extends passed the desired location. 
%We then obtain the unit vector of v and multiply it by the distance we want to get the new location

Q1=[cx,cy];

if isempty(Q2temp)
    scale_f=2;
    if dir==1
        Q2=[cx+(obs_length/2)*scale_f,cy+(obs_length/2)*scale_f]; %down to the right. length/2 ensures that we'll definitely be past the desired coordinate
    elseif dir==2
        Q2=[cx-(obs_length/2)*scale_f,cy+(obs_length/2)*scale_f]; %down to the left.
    end
else
    Q2=Q2temp;
end

V=Q2-Q1;
U=V/norm(V);
shift_dist=(obs_length/2)-P;

if T==1 %towards Q2
    R=Q1+(shift_dist*U);
elseif T==2 %away from Q2
    R=Q1-(shift_dist*U);
end

new_cx=R(1);
new_cy=R(2);

return



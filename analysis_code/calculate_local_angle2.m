close all;
%clearvars -except dat_all info_all
home

%Here we will claculate local changes in angles to determine how long the
%pEC should be for DTTs

%we will try this for different windows (ie samples to skip in the trajectory)

%do it for individual subjects, then pool the data together

win = 3;
% screen_offset = 0.1;
screen_offset = 0;

dtx = xbln_dt;
dty = ybln_dt;

ctgt = tgt(idx_all.dt_bln,:);
ctgt = ctgt(:,12); %which is left or right (0 is left)

tsx=24638/1000; %start target location for the tablet
tsy=27813/1000;
D=20;
left_30_tgt=[tsx+D*sind(-30),tsy+D*cosd(-30)]; %tablet target locations
right_30_tgt=[tsx+D*sind(30),tsy+D*cosd(30)];

[~, nt] = size(xbln_dt);

%OK, now loop through it
max_samples = 350;
local_ang = nan(num_subjects, nt, floor(max_samples / win));

ypos_max_angle = nan(num_subjects, nt);

for q=1:num_subjects
   
    for p=1:nt %looping through trials
        
        x_data = (dtx{q,p});
        y_data = dty{q,p};
        c=1;
         [ns] = length(x_data); 
         
         %for the current trial, know where the target is
         if ctgt(p)==0, tgt_loc = left_30_tgt; else, tgt_loc = right_30_tgt; end
         
        if sum(isnan(x_data)) ~= length(x_data)
            for k=win+1:win:ns-win
                %if p==2 & q==1, figure; hold on; plot(x_data, y_data); keyboard; end
                
                prev_x = x_data(k-win); %previous sample
                prev_y = y_data(k-win);
                
                cx = x_data(k); %current sample
                cy = y_data(k);
                
                nx = x_data(k+win);
                ny = y_data(k+win);
                
                v1 = [cx - prev_x, cy - prev_y];
                v2 = [nx - cx, ny - cy];
             
                angle = atan2d(det([v1;v2]),dot(v1,v2));
                
                %if angle< 1e-3, angle = nan; end
                
                local_ang(q,p,c) = abs(angle);
                c=c+1;
                 
                
            end
            
            %now lets do work on the angles for the current trial
            
            %turns out the angle calculation is a bit unstable at the
            %beginning of the movement where its smooth and close to 0
            %we will exclude data from the beginning and end of the trial
            %since its impossible for these to be relevant (the cue only
            %happens 2 cm into the movement)
            
            %if p==2 & q==2, keyboard; end
            
            displ = sqrt (x_data.^2 + y_data.^2);
            displ = displ - displ(1);
            idx1 = floor( find(displ <= 2, 1,'last') / win);
            %idx2 = floor( find(displ >= 18, 1,'first') / win); 
            %not the best way to handle this, 
            %excluding feedback corrections should be based on distance from target
            
            dist_tgt = sqrt( ((x_data) - ones(size(x_data)) * tgt_loc(1)).^2 +(y_data - ones(size(y_data)) * tgt_loc(2)).^2);
            idx2 = floor( find(dist_tgt <= 2,1,'first')/win);
            
            qq = squeeze(local_ang(q,p,:));
            
            qq(1:idx1) = nan; 
            %qq(idx2:end) = nan;
            
            [~, tmp1] = max(qq);
            tmp1 = tmp1 * win;
            
%             ypos_max_angle(q, p) = y_data(tmp1) - screen_offset;
              ypos_max_angle(q, p) = displ(tmp1);
            
        end
        
    end
end


yy = reshape(ypos_max_angle,nt*num_subjects,1);


hist(yy,50);
xlim([0,20]);
ylabel('number of occurrences');
xlabel('divergence point in y-position (cm)');

yd_sub = nanmedian(ypos_max_angle,2);
yd_sub_err = nanstd(ypos_max_angle,0,2);

figure; hold on;
plot([1:26],yd_sub,'bo');
display_xy_error_V2([1:26],yd_sub, [], yd_sub_err, 'k');

xlabel('subject');
ylabel('divergence point in y (cm)');
xlim([0,num_subjects+1]);
ylim([0,20]);

title('baseline 2TT');

P = 60;

%L = prctile(yy,P)

L1 = prctile(ypos_max_angle,P, 2)

L2 = nanmean(prctile(ypos_max_angle,P, 2))



% 
% if p==6,
%     figure; hold on; plot(x_data, y_data);
%     qq = squeeze(local_ang(q,p,:));
%     figure; plot(qq,'o');
%     %figure; title('diff'); plot(diff(qq));
%     [~,idx] = max(qq)
%     
%     keyboard;
%     
% end



% for qq = 1:num_subjects
%     for kn = 1:num_trials
%         
%         current_pos = Pos_sub_all{qq,kn};
%         initial_idx = find(current_pos >= initial_marker,1,'first');
%         %now use this to get the vector for angle calculations
%         x_data = xTraj_sub_all{qq,kn};
%         y_data = yTraj_sub_all{qq,kn};
%         
%         v1 = [x_data(1) - S(1),initial_marker]; %vector we will calculate the angle wrt
%         v2 = [x_data(initial_idx) - S(1), y_data(initial_idx) - S(2)]; %vector corresponding to movement
%         if ~isempty(v2)
%             angle = atan2d(det([v1;v2]),dot(v1,v2)); %CCCW is positive, so multiply by -1 (??)
%         else
% %             keyboard; %these are bad movements
%             angle= nan;
%         end
%         
%         initial_angle(qq,kn) = angle;
%        
%         initial_mov{qq,kn}(:,1) = x_data(1:initial_idx);
%         initial_mov{qq,kn}(:,2) = y_data(1:initial_idx);
%                 
%     end
% end




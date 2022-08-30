function P = determine_path_type_821(xTraj,yTraj,ST_idx,obs_idx,isLeft,isRight, is30,is45)
%this function takes in a given subject's trajectories and determines if the path used on
%each trial is the optimal one or a motor averaging one based on the obstacle location
%P is a vector of 0's (motor averaging indicator) or 1's (optimal)
%also need to know if its a leftward or rightward trial
%and the direction (30 vs 45 degrees)

%NOTE: this code should probably exclude trials where subjects actually hit the obstacle

[L,num_trials] = size(xTraj);
P=nan(1,num_trials);

%% construct the obstacle locations (copy-pasted from group analysis)
D=20; %20 cm movements
tsx=24638/1000; %start target location for the tablet
tsy=27813/1000;
S=[tsx,tsy];

left_30_tgt=[tsx+D*sind(-30),tsy+D*cosd(-30)]; %tablet target locations
right_30_tgt=[tsx+D*sind(30),tsy+D*cosd(30)];
left_45_tgt=[tsx+D*sind(-45),tsy+D*cosd(-45)];
right_45_tgt=[tsx+D*sind(45),tsy+D*cosd(45)];

%construct the obstacle centers for all cases
LT30_cent=((left_30_tgt-S)/2) + S;      %LT30x = LT30_cent(1);
RT30_cent=((right_30_tgt-S)/2) + S;     %RT30x = RT30_cent(1);
LT45_cent=((left_45_tgt-S)/2) + S;      %LT45x = LT45_cent(1);
RT45_cent=((right_45_tgt-S)/2) + S;     %RT45x = RT45_cent(1);

%%
for k=1:num_trials
    if ST_idx(k)==1 && obs_idx(k)==1 %only perform this analysis if we have a STT w/ an obstacle
        current_trial=[xTraj{:,k},yTraj{:,k}];
        x=xTraj{:,k};   y=yTraj{:,k};
        
        %determine which obstacle we're dealing with
        if isLeft(k)==1
            if is30(k)==1
                obsx = LT30_cent(1);
                obsy = LT30_cent(2);
            elseif is45(k)==1
                obsx = LT45_cent(1);
                obsy = LT45_cent(2);
            end
        elseif isRight(k)==1
            if is30(k)==1
                obsx = RT30_cent(1);
                obsy = RT30_cent(2);
            elseif is45(k)==1
                obsx = RT45_cent(1);
                obsy = RT45_cent(2);
            end
            
        else %straight ahead case
            obsx = nan;
        end
        
        if ~isnan(obsx)
            obs_loc=[obsx,obsy];
            %now find the closest point to compare
            distances = sqrt(sum(bsxfun(@minus, current_trial, obs_loc).^2,2));
            [~,closest_point_idx]= min(distances);
            closest_point=[x(closest_point_idx),y(closest_point_idx)];
            
            %use a reference vector to calculated singed angle?
            %         angle = atan2(closest_point(2), closest_point(1)) - atan2(obs_loc(2),obs_loc(1));
            %         if (angle < 0), angle=angle+  2 * pi; end
            %         angle=rad2deg(angle);
            %
            %         angle = atan2d(closest_point(2), 33) - atan2(obs_loc(2),obs_loc(1));
            
            if closest_point(1) > obsx & isRight(k)==1
                P(k) = 1;
            elseif closest_point(1) < obsx & isLeft(k)==1
                P(k) = 1;
            else
                P(k)=0;
            end
            
        end
        
    end
    
    
end





end
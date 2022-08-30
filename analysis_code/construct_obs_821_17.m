function [obs_rect]=construct_obs_821_17(obs_type,cx,cy,angle,tgt_dir, obs_length)

obs_centx=cx;
obs_centy=cy;

obs_width=1; %all obstacles have a width of 1 cm

if obs_type==1
    %obs_length=12;
    T=1;
    P=2;
elseif obs_type==2
    %obs_length=2;
    T=2;
    P=0;
elseif obs_type==4
    %obs_length=1;
    T=2;
    P=-obs_length;
    
elseif obs_type==99
    T = 1;
    P = 0;
    
end

[obs_temp,~,~]=make_obs_loc(obs_centx,obs_centy,angle,obs_width,obs_length);
edge_point=(obs_temp(1,:) + obs_temp(2,:))/2;
[ncx,ncy]=calculate_obstacle_offset(obs_centx,obs_centy,edge_point,P,obs_length,tgt_dir,T);
[P2,~,~]=make_obs_loc(ncx,ncy,angle,obs_width,obs_length);
obs_rect=[P2;P2(1,:)];
end

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

end
function [indp2, dep2] = transform_traj(indp1, dep1)
%%%in this function we will transform x and y trajectories into x(y) or y(x),
%and optionally, then transform back in terms of time

%%%this will take care of averaging trajectories together that have vastly
%different movement speeds

%%%We will do this by interpolating the dependent variable's trajectories
%%%onto a vector of the indepdent variable's positions every 0.254 mm
%%%(resolution of the tablet)





%first interpolate the dependent variable
max_samples = 300;
tmp1 = [0:0.0254: 0.0254 * (max_samples-1)];
%tmp1 = [min(indp1) : .0254 : max(indp1)];

indp2 = interp1([1:size(indp1,1)],indp1, linspace(indp1(1));

%interpolate the independent wrt the depedent
dep2 = interp1([1:size(dep1,1)],dep1, indp2);









end






% [tmpy, idx_y] = unique(y);
% tmpx = x(idx_y);
% plot(tmpx,tmpy)
% hold on; plot(x,y)
% qq = interp1(tmpy, tmpx, yy)
% plot(x,y,'--'); hold on; plot(qq,yy)



% num_trials=length(ydata);
% max_samples=max(cellfun(@(x) length(x), ydata));
% interp_data=zeros(max_samples,num_trials);

% %first interpolate all trials...
% for k=1:num_trials
%     current_trialx=cell2mat(xdata(k));    
%     current_trialy=cell2mat(ydata(k));
%     if sum(isnan(current_trialx))==length(current_trialx) %for some reason the tablet time vectors are nans...need to investigate
%         current_trialx=[0:0.005:0.005*(length(current_trialy)-1)];
%     end 
%     x_interp=linspace(current_trialx(1),current_trialx(end),max_samples);
%     interp_data(:,k)=interp1(current_trialx,current_trialy,x_interp);    
% end



max_samples = 174;
current_trialx = Time_sub_all{7,214};
current_trialy = xTraj_sub_all{7, 214};
x_interp=linspace(current_trialx(1),current_trialx(end),max_samples);

idx1 = interp1(current_trialx,current_trialy,x_interp);   



current_trialx = Time_sub_all{7,214};
current_trialy = yTraj_sub_all{7, 214};
x_interp=linspace(current_trialx(1),current_trialx(end),max_samples);

idy1 = interp1(current_trialx,current_trialy,x_interp);



xx = ( idx1' + xTraj_sub_all{7, 126} ) / 2;
yy = ( idy1' + yTraj_sub_all{7, 126} ) / 2;


plot(xx,yy,'.');


yy=[min(y) : .0254 : max(y)];
[tmpy, idx_x] = unique(x);
tmpx = y(idx_x);
plot(tmpx,tmpy)
%hold on; plot(x,y)
qq = interp1(tmpy, tmpx, yy)
plot(x,y,'--'); hold on; plot(qq,yy)
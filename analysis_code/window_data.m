function [new_data,std_mat]=window_data(xdata,ydata,win)
%interpolates and then windows data with a trial window of win
%x and y should have the same dimensions
%last window need not be a divisor

%NOTE: add condition for win=1

num_trials=length(ydata);
max_samples=max(cellfun(@(x) length(x), ydata));
interp_data=zeros(max_samples,num_trials);

%first interpolate all trials...
for k=1:num_trials
    current_trialx=cell2mat(xdata(k));    
    current_trialy=cell2mat(ydata(k));
    if sum(isnan(current_trialx))==length(current_trialx) %for some reason the tablet time vectors are nans...need to investigate
        current_trialx=[0:0.005:0.005*(length(current_trialy)-1)];
    end 
    x_interp=linspace(current_trialx(1),current_trialx(end),max_samples);
    interp_data(:,k)=interp1(current_trialx,current_trialy,x_interp);    
end
    
%now bin the data
if win~=1
    num_rem=rem(num_trials,win);
    if num_rem>0,
        win_flag=1;
        num_div=floor(num_trials/win);
        new_data=zeros(max_samples,num_div+1);
        std_mat=zeros(max_samples,num_div+1);
        for p=0:num_div-1,
            try,
                new_data(:,p+1)=nanmean(interp_data(:,p*win+1:(p+1)*win),2);
                std_mat(:,p+1)=nanstd(interp_data(:,p*win+1:(p+1)*win),0,2);
            catch, keyboard; end
        end
        new_data(:,end)=nanmean(interp_data(:,num_trials-num_rem+1:num_trials),2);
        std_mat(:,end)=nanstd(interp_data(:,num_trials-num_rem+1:num_trials),0,2);
    else
        win_flag=0;
        num_div=(num_trials/win);
        new_data=zeros(max_samples,num_div);
        std_mat=zeros(max_samples,num_div);
        for p=0:num_div-1,
            new_data(:,p+1)=nanmean(interp_data(:,p*win+1:(p+1)*win),2);
            std_mat(:,p+1)=nanstd(interp_data(:,p*win+1:(p+1)*win),0,2);
        end
    end
else
    new_data=interp_data;
    std_mat=nanstd(interp_data,0,2); 
end
return


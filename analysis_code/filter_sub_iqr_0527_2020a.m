function [filtered_data, bad_idx] = filter_sub_iqr_0527_2020a (data, iqr_num)

%the data will come in as vectors (for each subject)

[num_trials] = length(data);

bad_idx =zeros(1,num_trials); %let 1 represent a bad trial
filtered_data = data;

mean_angle = nanmedian(data);
current_iqr = iqr(data);

for p=1:num_trials
    current_sample =  data(p);
    
    if ~isnan(current_sample)
        
        if ( current_sample > mean_angle + (iqr_num * current_iqr) )  | ...
                ( current_sample < mean_angle - (iqr_num * current_iqr) )
            
            bad_idx(p) = 1;
            filtered_data(p) = nan;
        end
    end
end

end
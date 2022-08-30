function [filtered_data, bad_idx] = filter_iqr_0415a (data, iqr_num)

%the data will come in as a matrix for all subjects




[num_subjects, num_trials] = size(data);

bad_idx =zeros(num_subjects,num_trials); %let 1 represent a bad trial
filtered_data = data;

for q=1:num_subjects
    
    %d = data(q,:);
    mean_angle = nanmean(data(q,:));
    current_iqr = iqr(data(q,:));
    
    for p=1:num_trials        
        current_sample =  data(q,p);
        
        if ~isnan(current_sample)
            
            if ( current_sample > mean_angle + (iqr_num * current_iqr) )  | ...
                    ( current_sample < mean_angle - (iqr_num * current_iqr) )
                
                bad_idx(q,p) = 1;
                filtered_data(q,p) = nan;
            end
        end
    end
    
end

end
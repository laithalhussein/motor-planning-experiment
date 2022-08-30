function run_parse_0319_2018(Cdr)
close all;
delta_t = 0;

try
    if Cdr(end)~='\', Cdr=strcat(Cdr,'\'); end
catch
    keyboard;
end

 subject_list = {'aag', 'aiu', 'anm', 'asl', 'cah', 'ckw', 'cmg', 'dac' 'hns','iga', 'jcs', ...
     'lmk', 'lnc', 'maa', 'mav', 'mba','meo', 'mml', 'rab', 'rks', 'rlp', 'shr', 'tak', 'tdh', ...
     'tje', 'wme'};
[file_dates,num_sets] = get_dates(subject_list,Cdr);

DAT = cell(1,length(subject_list));

for n = 1:length(subject_list)
    tmp_dat = cell(1,num_sets(n));
    for s = 1:num_sets(n)
        f_name = [Cdr,subject_list{n},'_set_','a'+s-1,'_',file_dates{n},'.mat']; %CHANGE LETTER ACCORDINGLY
        tmp_dat{s} = load(f_name);
    end
    DAT{n} = parse_VMR_subject_obs_V5(tmp_dat,delta_t);
end

% if n==1, save_name=['S_',subject_list{n},'_',file_dates{n},'.mat'];
% else, save_name=['obs_exp_dat_all_',date,'.mat']; end

save_name=['obs_exp_parse_all_',date,'.mat'];

save(save_name,'DAT');

end


function [file_dates,num_sets] = get_dates(f_names,dir_)

for i = 1:length(f_names)
    tmp = dir([dir_,f_names{i},'_set*']);
    num_sets(i) = length(tmp);
    try
    file_dates{i} = tmp(1).name(length(tmp(1).name)-14:length(tmp(1).name)-4);
    catch
        keyboard;
    end
end

end



load('st_bln_data');              qoDat=st_bln_data'; qo = stats2(qoDat);    % st_bln_data = expt 2b no-obstacle data: 6x8  cell all 1TTs (Nsubjects x [Left_test_epoch, Right, Center, Left_baseline, Right, Center])' ])'
load('st_bln_data_exp2a'); qaoDat=st_bln_data_exp2a'; qao = stats2(qaoDat);  % st_bln_data_exp2a = expt 2b no-obstacle data: 6x26 cell all 1TTs (Nsubjects x [Left_test_epoch, Right, Center, Left_baseline, Right, Center])' ])'
load('sm_data');         qDat=sm_data';         q = stats2(qDat);          % sm_data = expt 2b obstacle data:  4x26 cell (Nsubjects x [1TT-leftObstacle,1TT-right, 2TT-left, 2TT-right])'
load('exp2a_sm_data');   qaDat=exp2a_sm_data';  qa = stats2(qaDat); 


dtL = cell2mat(qDat(:,3)); %subject x trial
dtR = cell2mat(qDat(:,4));

stL = cell2mat(qDat(:,1));
stR = cell2mat(qDat(:,2));

Nb = size(dtL,1);
subs = [1:Nb];

%% plot mean vs median for each subject
close all;
grey = ones(3,1)/2;
figure;
subplot(221); hold on;
title('Left 2TT data');
dtL_mean = nanmean(dtL,2);
dtL_median = nanmedian(dtL,2);
plot(subs,dtL_mean,'color','k','Marker','o','Linestyle','None', 'displayname','Mean');
plot(subs,dtL_median,'color','r','Marker','o','Linestyle','None', 'displayname','Median');
xlabel('Subjects');
ylabel('Mov direction (deg)');
legend('show');

subplot(222); hold on;
title('Right 2TT data');
dtR_mean = nanmean(dtR,2);
dtR_median = nanmedian(dtR,2);
plot(subs,dtR_mean,'color','k','Marker','o','Linestyle','None', 'displayname','Mean');
plot(subs,dtR_median,'color','r','Marker','o','Linestyle','None', 'displayname','Median');
xlabel('Subjects');
ylabel('Mov direction (deg)');

subplot(223); hold on;
title('Left 1TT data');
stL_mean = nanmean(stL,2);
stL_median = nanmedian(stL,2);
plot(subs,stL_mean,'color','k','Marker','o','Linestyle','None', 'displayname','Mean');
plot(subs,stL_median,'color','r','Marker','o','Linestyle','None', 'displayname','Median');
xlabel('Subjects');
ylabel('Mov direction (deg)');

subplot(224); hold on;
title('Right 1TT data');
stR_mean = nanmean(stR,2);
stR_median = nanmedian(stR,2);
plot(subs,stR_mean,'color','k','Marker','o','Linestyle','None', 'displayname','Mean');
plot(subs,stR_median,'color','r','Marker','o','Linestyle','None', 'displayname','Median');
xlabel('Subjects');
ylabel('Mov direction (deg)');

%determine difference between mean and median
dtL_diff = abs(dtL_mean-dtL_median);
dtR_diff = abs(dtR_mean-dtR_median);
stL_diff = abs(stL_mean-stL_median);
stR_diff = abs(stR_mean-stR_median);

figure; hold on;
title('Difference between mean and median');
plot(subs,dtL_diff, 'k','displayname', 'left 2TT');
plot(subs,dtR_diff, 'color',grey,'displayname', 'right 2TT');
plot(subs,stL_diff, 'b','displayname', 'left 1TT');
plot(subs,stR_diff, 'r','displayname', 'right 1TT');
xlabel('subjects');
ylabel('absolute difference between mean and median (deg)');
legend('show');

%data appear pretty well behaved, except for 1 subject in the right 2TT
%(21) and 2 subjects in left 1TT (16 & 25) -- lets investigate

figure; hold on;
title('right 2TT: subject 21');
plot(dtR(21,:)); %3-5 trials at most that seem like outliers??

figure; hold on;
title('left 1TT: subject 16');
plot(stL(16,:)); %4 trials at most that seem like obvious outliers???

figure; hold on;
title('left 1TT: subject 25');
plot(stL(25,:)); %2 trials at most that seem potential outliers

%%%check to see if data fall within IQR threshold
iqr_thresh = 2.5;
% bad_trial = zeros(

%%%plot subject data used in assymettry analysis


% plot(nanmedian(dtL,2),nanmean(dtL,2),'color',grey,'Marker','o','Linestyle','None');
% plot([-5:15],[-5:15],'k--');
% axis tight;




%%
function y = stats2(q)
for k1=1:size(q,1), 
    for k2=1:size(q,2), 
        d=q{k1,k2}'; 
        y.N(k1,k2)=length(d);  y.N1(k1,k2)=sum(isfinite(d)); 
        y.M(k1,k2)=nanmean(d); y.Md(k1,k2)=nanmedian(d); 
        y.S(k1,k2)=nanstd(d);  y.Sd(k1,k2)=nanstd(dtrend(d,11)); 
        y.I(k1,k2)=iqr(d);     y.Id(k1,k2)=iqr(dtrend(d,11));
    end
    if k2==4,
        qt{k1,1} = [q{k1,1}, q{k1,2}];
        qt{k1,2} = [q{k1,3}, q{k1,4}];
        qt1{k1,1} = [q{k1,1}-nanmean(q{k1,1}), q{k1,2}-nanmean(q{k1,2})];
        qt1{k1,2} = [q{k1,3}-nanmean(q{k1,3}), q{k1,4}-nanmean(q{k1,4})];
    elseif k2==6,
        qt{k1,1} = [q{k1,1}, q{k1,2}];
        qt{k1,2} = [q{k1,4}, q{k1,5}];
        qt1{k1,1} = [q{k1,1}-nanmean(q{k1,1}), q{k1,2}-nanmean(q{k1,2})];
        qt1{k1,2} = [q{k1,4}-nanmean(q{k1,4}), q{k1,5}-nanmean(q{k1,5})];
    end
    
    
    for k2=1:size(qt,2),
        d=qt{k1,k2}'; 
        y.N1(k1,k2)=length(d);  y.N11(k1,k2)=sum(isfinite(d)); 
        y.M1(k1,k2)=nanmean(d); y.Md1(k1,k2)=nanmedian(d); 
        d1=qt1{k1,k2}'; 
        y.S1(k1,k2)=nanstd(d1);  y.Sd1(k1,k2)=nanstd(dtrend(d1,11)); 
        y.I1(k1,k2)=iqr(d1);     y.Id1(k1,k2)=iqr(dtrend(d1,11)); 
    end
end

function y = dtrend(x,n)
for k=1:size(x,2), y(:,k) = x(:,k)-smooth(x(:,k),n); end
end
end


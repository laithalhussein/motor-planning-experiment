qq = (-s_AL + s_AR)/2;

b1 = nan(26,1);
f = 0;
for k=1:26
    for kk = 1:46
        if abs(qq(k,kk))/2>(30-phiA)
            f = f+1;
        end
    end
    
    b1(k) = f;
    f = 0;
end


figure; hold on;

plot([1:26], b1/46*100, '.-', 'markersize', 10, 'color', 'k');

xlabel('subject #');
ylabel('fraction of predictions');
title('fraction of MA predictions that exceed the obstacle bound based on single trial estimates of 1-target trial IMD');
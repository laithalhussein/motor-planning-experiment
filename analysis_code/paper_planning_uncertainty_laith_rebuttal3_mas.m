load('st_bln_data');     
qoDat=st_bln_data';    
qo = stats2(qoDat);

load('sm_data');         
qDat=sm_data';         
q = stats2(qDat);

load('exp2a_sm_data');   
qaDat=exp2a_sm_data';  
qa = stats2(qaDat);

%%

[e2b_r, e2b_p] = corr(q.S1(:,2),q.M1(:,2));
[e2a_r, e2a_p] = corr(qa.S1(:,2),qa.M1(:,2));

[~,p]=corr([qo.I(:,4:6),q.I],'rows','pairwise'); round(-log10(p))
%[~,p]=corr([qo.Id(:,4:6),q.Id],'rows','pairwise'); round(-log10(p))
%disp(['Correlation between SM and MV for 2TTs in expt 2b is:  ', num2str()]);



% Fig R3:  Direction-specificity within 1TTs expt 2b

figure;  % Fig R3:  Direction-specificity of the relationship between SM & VAR for left vs right 1TTs in expt2b (~100deg separation for L vs R)
subplot(131); [rc,pc,pd] = plotxy1([q.S(:,1);q.S(:,2)], -[q.M(:,2);q.M(:,1)]);  % negative sign due to orientation convention used
title({'Modulation of SM by other-side variability',['Left vs right 1-target trials:   r = ',num2strFP(rc,2), ',  p = ',num2strFP(pc,2)]})
xlabel('Observed other-side variability');
ylabel('Observed safety margin'); grid
legend({'indivudal side data','group mean','group median'})

subplot(132); [rc,pc,pd] = plotxy1([q.S(:,1);q.S(:,2)], -[q.M(:,1);q.M(:,2)]);  % negative sign due to orientation convention used
title({'Evidence for intra-individual modulation of SM by variability',['Left vs right 1-target trials:   r = ',num2strFP(rc,2), ',  p = ',num2strFP(pc,2)]})
xlabel('Observed same-side variability');
ylabel('Observed safety margin'); grid
legend({'indivudal side data','group mean','group median'})

subplot(133); [rc,pc,pd] = plotxy1(100*(q.S(:,1)-q.S(:,2))./q.S(:,2),100*(q.M(:,1)-q.M(:,2))./q.M(:,2)); 
title({'Evidence for intra-individual modulation of SM by variability',['Left vs right 1-target trials:   r = ',num2strFP(rc,2), ',  p = ',num2strFP(pc,2)]})
xlabel('Percent change in variability from right to left');
ylabel('Percent change in safety margin from right to left'); grid
legend({'indivudal particpant data','group mean','group median'})
set(gcf,'position', [0 0 1555 404])


% % Direction-specificity within 1TTs expt 2a
% 
% figure;  % Fig ??:  Direction-specificity of the relationship between SM & VAR for left vs right 1TTs in expt2a (~20deg separation for L vs R)
% subplot(131); [rc,pc,pd] = plotxy1([qa.S(:,1);qa.S(:,2)], -[qa.M(:,2);qa.M(:,1)])  % negative sign due to orientation convention used
% title({'Modulation of SM by other-side variability',['Left vs right 1-target trials:   r = ',num2str(rc), ',  p = ',num2str(pc)]})
% xlabel('Observed other-side variability');
% ylabel('Observed safety margin'); grid
% legend({'indivudal side data','group mean','group median'})
% 
% subplot(132); [rc,pc,pd] = plotxy1([qa.S(:,1);qa.S(:,2)], -[qa.M(:,1);qa.M(:,2)])  % negative sign due to orientation convention used
% title({'Evidence for intra-individual modulation of SM by variability',['Left vs right 1-target trials:   r = ',num2str(rc), ',  p = ',num2str(pc)]})
% xlabel('Observed same-side variability');
% ylabel('Observed safety margin'); grid
% legend({'indivudal side data','group mean','group median'})
% 
% subplot(133); [rc,pc,pd] = plotxy1(100*(qa.S(:,1)-qa.S(:,2))./qa.S(:,2),100*(qa.M(:,1)-qa.M(:,2))./qa.M(:,2));
% title({'Evidence for intra-individual modulation of SM by variability',['Left vs right 1-target trials:   r = ',num2str(rc), ',  p = ',num2str(pc)]})
% xlabel('Percent change in variability from right to left');
% ylabel('Percent change in safety margin from right to left'); grid
% legend({'indivudal particpant data','group mean','group median'})



% Fig R4:  Direction-specificity 2TTs vs 1TTs expt 2b

figure;   % Fig R4:  Direction-specificity of the relationship between SM & VAR for 2TTs vs 1TTs in expt2b 
subplot(221); [rc,pc,pd] = plotxy1(q.S1(:,1), -q.M1(:,1))
title({'Relationship between SM and MV for 1TTs [spatially matched]',['corr = ',num2strFP(rc,2), ',  p_c_o_r_r = ',num2strFP(pc,2)]});
xlabel('1-target trial (1TT) variability (SD) [deg]');
ylabel('1-target trial (1TT) safety margin [deg]'); grid
legend({'individual particpant data','group mean','group median'})

subplot(222); [rc,pc,pd] = plotxy1(q.S1(:,2), -q.M1(:,1))
title({'Relationship between 1TT SM and 2TT MV [spatially mismatched]',['corr = ',num2strFP(rc,2), ',  p_c_o_r_r = ',num2strFP(pc,2)]});
xlabel('2-target trial (2TT) variability (SD) [deg]');
ylabel('1-target trial (1TT) safety margin [deg]'); grid
legend({'individual particpant data','group mean','group median'})

subplot(223); [rc,pc,pd] = plotxy1(q.S1(:,1),15+q.M1(:,2));
title({'Relationship between 2TT SM and 1TT MV [spatially mismatched]',['corr = ',num2strFP(rc,2), ',  p_c_o_r_r = ',num2strFP(pc,2)]});
xlabel('1-target trial (1TT) variability (SD) [deg]');
ylabel('2-target trial (2TT) safety margin [deg]'); grid
legend({'individual particpant data','group mean','group median'})

subplot(224); [rc,pc,pd] = plotxy1(q.S1(:,2),15+q.M1(:,2));
title({'Relationship between SM size and MV for 2TTs [spatially matched]',['corr = ',num2strFP(rc,2), ',  p_c_o_r_r = ',num2strFP(pc,2)]});
xlabel('2-target trial (2TT) variability (SD) [deg]');
ylabel('2-target trial (2TT) safety margin [deg]'); grid
legend({'individual particpant data','group mean','group median'})
set(gcf,'position', [0   0   850   732])



% Fig R5:   Direction-specificity of the relationship between delta SM & delta VAR between 2TTs & 1TTs in expt2b 

figure; % Fig R4:  Direction-specificity of the relationship between SM & VAR for 2TTs vs 1TTs in expt2b 
subplot(131); [rc,pc,pd] = plotxy1((q.S(:,1)-q.S1(:,2))./q.S1(:,2), (-q.M(:,1)-(15+q.M1(:,2)))./(15+q.M1(:,2)),[]);
title({'Comparison between 2TTs and left-side 1TTs',['corr = ',num2strFP(rc,2), ',  p_c_o_r_r = ',num2strFP(pc,2)]})
xlabel('% Change in left-side 1-target trial (1TT) variability');
ylabel('% Change in 2-target trial (2TT) safety margin (SM)');
grid
legend({'individual particpant data','group mean','group median'})

subplot(132); [rc,pc,pd] = plotxy1((q.S(:,2)-q.S1(:,2))./q.S1(:,2), (-q.M(:,2)-(15+q.M1(:,2)))./(15+q.M1(:,2)),[]);
title({'Comparison between 2TTs and right-side 1TTs',['corr = ',num2strFP(rc,2), ',  p_c_o_r_r = ',num2strFP(pc,2)]})
xlabel('% Change in right-side 1-target trial (1TT) variability');
ylabel('% Change in 2-target trial (2TT) safety margin (SM)');
grid
legend({'individual particpant data','group mean','group median'})

subplot(133); [rc,pc,pd] = plotxy1((q.S1(:,1)-q.S1(:,2))./q.S1(:,2), (-q.M1(:,1)-(15+q.M1(:,2)))./(15+q.M1(:,2)),[]);
title({'Comparison between 2TTs and combined 1TT data',['corr = ',num2strFP(rc,2), ',  p_c_o_r_r = ',num2strFP(pc,2)]})
xlabel('% Change in pooled 1-target trial (1TT) variability');
ylabel('% Change in 2-target trial (2TT) safety margin (SM)');
grid
legend({'individual particpant data','group mean','group median'})
set(gcf,'position', [0 0 1555 404])



% Fig R6:  Direction-specificity for wide vs narrow spatial separation (comparing expt 2b vs 2a)

figure; % Left vs Right within the 1TT data
a=q.S(:,1:2); a=a(:); b=q.M(:,1:2); b=b(:);  b2=q.M(:,[2 1]); b2=b2(:); c1 = [corr(a,-b); corr(a,-b2)];
a=qa.S(:,1:2); a=a(:); b=qa.M(:,1:2); b=b(:);  b2=qa.M(:,[2 1]); b2=b2(:); ca1 = [corr(a,-b); corr(a,-b2)];
cc=[c1;-ca1];
green = [0 0.5 0];  gray = [1 1 1]*0.8;  
bar([1 4],cc([1 3]),0.25,'FaceColor', green); hold on;
bar([2 5],cc([2 4]),0.25,'FaceColor', gray);
%title({'Comparison showing that unmatched safety margin (SM) and motor variability (MV) leads to','greater reductions in prediction quality when spatial separation is wide, than when it is narrow'} );
%legend('Matched side','Unmatched side')
text(1.5, 0.71, {'Expt 2b data (from Fig R3)','(wide spatial separation)','Left vs Right 1TT data'},'horizontalAlignment','center','fontsize',12,'fontweight','bold')
text(4.5, 0.71, {'Exp1 2a data','(narrow spatial separation)','Left vs Right 1TT data'},'horizontalAlignment','center','fontsize',12,'fontweight','bold')
label1 = {'MV predicts same-side SM','MV predicts opposite-side SM'};
% xticks([1 2 4 5]);
% xticklabels({label1{:},label1{:}})
% xtickangle(45)
% ylabel('Correlation coefficient')
% ylim([-.15 0.89])
% analysis for 1TT vs 2TT data
c1 = corr(q.S1(:,1:2),q.M1(:,1:2).*[-1,1]); c1=c1(:); c1=c1([1 4 2 3]);
ca1 = corr(qa.S1(:,1:2),qa.M1(:,1:2)); ca1=ca1(:); ca1=ca1([1 4 2 3]);
cc=[c1;ca1];
green = [0 0.5 0];  gray = [1 1 1]*0.8;  
bar(7+[1 2 6 7],cc([1 2 5 6]),'FaceColor', green); hold on;
bar(7+[3 4 8 9],cc([3 4 7 8]),'FaceColor', gray);
title({'Comparison showing that unmatched safety margin (SM) and motor variability (MV) leads to','greater reductions in prediction quality when spatial separation is wide than when it is narrow, consistent with spatial secificity'} );
legend('Matched side','Unmatched side')
text(7+2.5, 0.81, {'Expt 2b data (from Fig R4)','(wide spatial separation)','1TT vs 2TT data'},'horizontalAlignment','center','fontsize',12,'fontweight','bold')
text(7+7.5, 0.81, {'Expt 2a data','(narrow spatial separation)','1TT vs 2TT data'},'horizontalAlignment','center','fontsize',12,'fontweight','bold')
label2 = {'1TT MV predicts 1TT SM','2TT MV predicts 2TT SM','1TT MV predicts 2TT SM','2TT MV predicts 1TT SM'};
xticks([1:4,6:9]);
xticks([1,2,4,5,7+[1:4,6:9]]);
xticklabels({label1{:},label1{:},label2{:},label2{:}})
xtickangle(45)
ylabel('Correlation coefficient')
ylim([-.18 0.92])
set(gcf,'position', [0 0 1286 441])

%% indivoduation analysis

po{1} = mean(-q.M1(:,1)).*mean(q.S1(:,2))./mean(q.S1(:,1)).*q.M1(:,1)./q.M1(:,1);
po{2} = (-q.M1(:,1)).*mean(q.S1(:,2))./mean(q.S1(:,1));
po{3} = mean(-q.M1(:,1)).*(q.S1(:,2))./mean(q.S1(:,1));
po{4} = mean(-q.M1(:,1)).*mean(q.S1(:,2))./(q.S1(:,1));
po{5} = (-q.M1(:,1)).*(q.S1(:,2))./(q.S1(:,1));
ma{1} = mean(q.M1(:,1))/2;
ma{2} = (q.M1(:,1))/2;

clear rmsPo rmsMa
for k1=1:100, for k=1:length(po), rmsPo(k1,k) = rms(q.M1(:,2)-k1/100*(po{k}-15)); end; end
for k1=1:100, for k=1:length(ma), rmsMa(k1,k) = rms(q.M1(:,2)-k1/100*(ma{k}-15)); end; end
clear sdPo sdMa
for k1=1:100, for k=1:length(po), sdPo(k1,k) = std(q.M1(:,2)-k1/100*(po{k}-15)); end; end
for k1=1:100, for k=1:length(ma), sdMa(k1,k) = std(q.M1(:,2)-k1/100*(ma{k}-15)); end; end
figure; plot(rmsPo)
figure; plot(rmsMa)
figure; plot(sdPo)
figure; plot(sdMa)
xlabel('Obstacle avoidance weighting (%)');
ylabel('RMS (thin) or SD (bold) of individual partiipsnt error (deg)')
legend({'Group Mean', 'Individual M1', 'Individual S2', 'Individual S1', 'Fully Individualized'})

%%

reg1 = [-q.M1(:,1),q.S1(:,1),q.S1(:,2)];
figure; for k=1:3, subplot(1,3,k); [rc,pc,pd] = plotxy1(q.M1(:,2),reg1(:,k)); title(num2str([rc,pc])); end
[B,BINT,R,RINT,STATS] = regress(q.M1(:,2),[reg1,reg1(:,1)*0+1]); STATS
[B,BINT,R,RINT,STATS] = regress(log(q.M1(:,2)),[reg1,reg1(:,1)*0+1]); STATS

keyboard

figure; 
for k=1:5, subplot(1,5,k); plotxy1((po{k}-15)/2,q.M1(:,2),1); title(  std(q.M1(:,2) - (po{k}-15)/2)  ); end

rmsPo
rmsMa

%%

% 
% keyboard
% 
% figure;   % Comparison of variability on 1TTs and 2TTs using both SD & IQR
% subplot(121); [rc,pc,pd] = plotxy1(q1S(:,1), q1S(:,2))
% title({'Data show no consistent increase in variability from 1-target to 2-target trials',['p = ',num2str(pd), '   corr = ',num2str(rc), ',  p_c_o_r_r = ',num2str(pc)]})
% xlabel('1-target trial variability (SD) [deg]');
% ylabel('2-target trial variability (SD) [deg]'); grid
% legend({'individual particpant data','group mean','group median','y=x'})
% 
% subplot(122); [rc,pc,pd] = plotxy1(q1I(:,1), q1I(:,2))
% title({'Data show no consistent increase in variability from 1-target to 2-target trials',['p = ',num2str(pd), '   corrr = ',num2str(rc), ',  p_c_o_r_r = ',num2str(pc)]})
% xlabel('1-target trial variability (IQR) [deg]');
% ylabel('2-target trial variability (IQR) [deg]'); grid
% legend({'individual particpant data','group mean','group median','y=x'})
% 
% 
% figure;  % Normalized changes in SM size vs normalized changes in variability (SD) for R vs L obstacle positions
% subplot(121); [rc,pc,pd] = plotxy1(qS(:,1)-qS(:,2), -(qM(:,1)-qM(:,2)))  % negative sign due to orientation convention used
% title({'Evidence for intra-individual modulation of SM by variability',['Left vs right 1-target trials:   r = ',num2str(rc), ',  p = ',num2str(pc)]})
% xlabel('Change in variability from right to left');
% ylabel('Change in safety margin from right to left'); grid
% legend({'indivudal particpant data','group mean','group median'})
% 
% subplot(122); [rc,pc,pd] = plotxy1(100*(qS(:,1)-qS(:,2))./qS(:,2),100*(qM(:,1)-qM(:,2))./qM(:,2));
% title({'Evidence for intra-individual modulation of SM by variability',['Left vs right 1-target trials:   r = ',num2str(rc), ',  p = ',num2str(pc)]})
% xlabel('Percent change in variability from right to left');
% ylabel('Percent change in safety margin from right to left'); grid
% legend({'indivudal particpant data','group mean','group median'})
% 
% 
% figure;   % Compare std computed using properly vs improperly pooled data
% subplot(121); [rc,pc,pd] = plotxy1(q1S(:,2), 15+qtM(:,2),1)
% title({'Using pooled variability from L/R data from 2TTs',['corr = ',num2str(rc), ',  p_c_o_r_r = ',num2str(pc)]})
% xlabel('2-target trial variability (SD) [deg]');
% ylabel('2-target trial safety margin [deg]'); grid
% legend({'individual particpant data','group mean','group median','y=x','y=mx+b'})
% 
% subplot(122); [rc,pc,pd] = plotxy1(qtS(:,2), 15+qtM(:,2),1);
% title({'Using pooled variability from L/R data from 2TTs combined WITHOUT taking means into account ',['corr = ',num2str(rc), ',  p_c_o_r_r = ',num2str(pc)]})
% xlabel('2-target trial variability (SD) [deg]');
% ylabel('2-target trial safety margin [deg]'); grid
% legend({'individual particpant data','group mean','group median','y=x','y=mx+b'})
% 
% 
% % Asymmetry analysis looking for decision variability re-routing
% figure;
% dS = sqrt(q1S(:,2).^2 - q1S(:,1).^2);
% dS1 = sort(dS(imag(dS)==0));
% a = real(dS).^2./q1S(:,1).^2;
% i = find(a > 0.5);
% [~,i1] = sort(a(i),'descend');
% N=1e6; Nh = [-35:1:20];
% for k=1:length(i), ii=i(i1(k));
%     xd0 = dS(ii)*randn(N,1)-qtM(ii,2);
%     MarginPoint = 15-2.87*q1S(ii,1);
%     xd = xd0; xd(xd>MarginPoint) = MarginPoint;
%     x = xd + q1S(ii,1)*randn(N,1);
%     x0 = xd0 + q1S(ii,1)*randn(N,1);
%     [hy0,hx0] = hist(x0,Nh);
%     [hy1,hx1] = hist(xd,Nh);
%     [hy2,hx2] = hist(x,Nh);
%     subplot(3,2,k); bar(hx1,hy1/1e6,'r'); hold on; plot(hx2,hy2/1e6,'b','linewidth',2); plot(hx0,hy0/1e6,'k','linewidth',2)
%     asym(1) = std(x0(x0>median(x0))) ./ std(x0(x0<median(x0)));
%     mxd=median(xd); asym(2) = (mean([xd(xd<mxd)-mxd].^2))^(1/2) ./ (mean([xd(xd>mxd)-mxd].^2))^(1/2); 
%     mxd=median(xd); asym(3) = (mean([xd(xd<mxd)-mxd].^4))^(1/4) ./ (mean([xd(xd>mxd)-mxd].^4))^(1/4); 
%     title (['decVar / motVar ratio:   ',num2str(real(dS(ii)).^2./q1S(ii,1).^2,3), '   Decision Asymmetry: ', num2str(asym(2),3), '   Output Asymmetry: ', num2str(asym(3),3)] )
%     ylim([0 0.135])
% end
% subplot(3,2,5); xlabel('Movement direction (deg)');
% subplot(3,2,6); xlabel('Movement direction (deg)');
% subplot(3,2,3); ylabel('Probabilty density');
% 
function s = num2strFP(x, n)  %fixed point (FP) improvement on num2str
extraZeros = -floor(log10(x))-1;
s = num2str(x,['%0.',num2str(n+extraZeros),'f']);
end

function y = stats1(q)
for k1=1:size(q,1), for k2=1:size(q,2), 
        d=q{k1,k2}'; 
        y.N(k1,k2)=length(d);  y.N1(k1,k2)=sum(isfinite(d)); 
        y.M(k1,k2)=nanmean(d); y.Md(k1,k2)=nanmedian(d); 
        y.S(k1,k2)=nanstd(d);  y.Sd(k1,k2)=nanstd(dtrend(d,11)); 
        y.I(k1,k2)=iqr(d);     y.Id(k1,k2)=iqr(dtrend(d,11)); 
    end; end
end

function y = stats2(q)
for k1=1:size(q,1), 
    for k2=1:size(q,2), 
        d=q{k1,k2}'; 
        y.N(k1,k2)=length(d);  y.N1(k1,k2)=sum(isfinite(d)); 
        y.M(k1,k2)=nanmean(d); y.Md(k1,k2)=nanmedian(d); 
        y.S(k1,k2)=nanstd(d);  y.Sd(k1,k2)=nanstd(dtrend(d,11)); 
        y.I(k1,k2)=iqr(d);     y.Id(k1,k2)=iqr(dtrend(d,11)); 
    end
    qt{k1,1} = [q{k1,1}, q{k1,2}];%concatinates left (flipped) and right 1TT
    qt{k1,2} = [q{k1,3}, q{k1,4}]; %concatinates left (flipped) and right 2TT
    qt1{k1,1} = [q{k1,1}-nanmean(q{k1,1}), q{k1,2}-nanmean(q{k1,2})]; 
    qt1{k1,2} = [q{k1,3}-nanmean(q{k1,3}), q{k1,4}-nanmean(q{k1,4})];
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


function [rc,pc,pd,B] = plotxy1(a,b,k)
if nargin < 3, k=[]; end

plot(a,b,'o'); hold on;
[rc,pc] = corr(a,b); 
[~,pd] = ttest(a,b);
plot(mean(a),mean(b),'sr','linewidth',2,'markersize',10);  
plot(median(a),median(b),'^r','linewidth',2,'markersize',10); 
a1=min([a(:)]);
a2=max([a(:)]);
a1 = a1 - 0.05*(a2-a1);
a2 = a2 + 0.05*(a2-a1);
if length(k)>0, plot([a1 a2],k*[a1 a2],'--k'); end
[B] = regress(b,[a,a*0+1]); plot([a1 a2],B(1)*[a1 a2]+B(2),'--b')
end

figure;
dS = sqrt(q1S(:,2).^2 - q1S(:,1).^2);
a = real(dS).^2./q1S(:,1).^2;
i = find(a > 0.5);
[~,i1] = sort(a(i),'descend');
N=1e6; Nh = [-35:1:20];
for k=1:length(i), ii=i(i1(k));
    xd0 = dS(ii)*randn(N,1)-qtM(ii,2);
    MarginPoint = 15-2.87*q1S(ii,1);
    xd = xd0; xd(xd>MarginPoint) = MarginPoint;
    x = xd + q1S(ii,1)*randn(N,1);
    x0 = xd0 + q1S(ii,1)*randn(N,1);
    [hy0,hx0] = hist(x0,Nh);
    [hy1,hx1] = hist(xd,Nh);
    [hy2,hx2] = hist(x,Nh);
    subplot(3,2,k); bar(hx1,hy1/1e6,'r'); hold on; 
    plot(hx2,hy2/1e6,'b','linewidth',2); plot(hx0,hy0/1e6,'k','linewidth',2)
    asym(1) = std(x0(x0>median(x0))) ./ std(x0(x0<median(x0)));
    asym(2) = std(xd(xd>median(xd))) ./ std(xd(xd<median(xd)));
    asym(3) = std(x(x>median(x))) ./ std(x(x<median(x)));
    title (['decVar / motVar ratio:   ',num2str(real(dS(ii)).^2./q1S(ii,1).^2,3), '   Decision Asymmetry: ', num2str(asym(2),3), '   Output Asymmetry: ', num2str(asym(3),3)] )
    ylim([0 0.135])
end
subplot(3,2,5); xlabel('Movement direction (deg)');
subplot(3,2,6); xlabel('Movement direction (deg)');
subplot(3,2,3); ylabel('Probabilty density');
function create_hist_init_sc(data,color,b)

% if isempty(b)
%     bins = 18;
% else 
%     bins = b;
% end

bins = 18;
alpha_val = 0.75;

if isempty(data), return, end

h = histogram(data,bins,'facecolor',color,'facealpha',alpha_val,'displaystyle','stairs');
h.EdgeColor = color;
h.LineWidth = 1;

%add vertical line @ mean
%mu = nanmean(data);
mu = nanmedian(data);
YL = ylim;

% h1=plot([mu, mu], [YL(1), YL(2)], 'color', ones(1,3)*0, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% leg1=legend(h1);
% set(leg1,'Location','Best');
%plot([T, T], [YL(1), YL(2)], 'color',ones(1,3) * 0, 'linestyle','-','linewidth',3);
% leg1 = legend([num2str(mu)]);
% set(leg1,'Location','Best');

%axis tight;
end
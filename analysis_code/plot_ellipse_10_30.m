function [ h, orient_deg ] = plot_ellipse_10_30(x,y,col,LW, error_type)

%need to plot an ellipse where the axes are defined by the noise we have in the data

num_errors = 1; %# of sd or se
num_samples = length(x(:));

c = nancov([x(:),y(:)]);
%get the eigenvalues/vectors
[evec,eval] = eig(c);

%start to define the ellipse points
a = [0:360]' * pi/180;

if regexpi(error_type, 'standard_deviation')
    e_points = [cos(a),sin(a)]*sqrt(eval)*evec' * num_errors;
    H=sqrt(eval)*evec' * num_errors;
elseif regexpi(error_type, 'standard_error')
    e_points = [cos(a),sin(a)]*sqrt(eval)*evec' * [num_errors / sqrt(num_samples)];
    H=sqrt(eval)*evec' * num_errors * [num_errors / sqrt(num_samples)];
end

h = plot(e_points(:,1)+nanmean(x), e_points(:,2)+nanmean(y),'color',col,'linewidth',LW); %

%h = plot(exy(:,1)+x_center, exy(:,2)+y_center);

uistack(h, 'bottom');
hg=get(h,'Annotation');
hLegendEntry = get(hg,'LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off');

% find the orientation
[max_eval,NDX]=max(max(abs(eval)));

orient_deg = 180/pi*atan2(evec(2,NDX),evec(1,NDX));

% get the area???
r_1=sqrt(H(1,1)^2+H(1,2)^2);
r_2=sqrt(H(2,1)^2+H(2,2)^2);
R=[r_1,r_2];

% ellipse_area = r_1*r_2*pi;

hold on;

%we could show the eigenvector itself...
%h1=plot([mean(x),mean(x)+R(NDX)*cosd(orient_deg)],[mean(y),mean(y)+R(NDX)*sind(orient_deg)],'LineWidth',3,'color','g');

shg
end
function p = determine_ma_pred(sm, bias, slope)

if isempty(slope), slope = 0.5; end

%find the MA prediction based on individual subject data and the population average

ma_gm =  slope*(-30 + sm) + (1-slope)*(30+mean(bias)) ;
ma_sub =  slope*(-30 + sm) + (1-slope)*(30+bias) ;

% ma_gm = slope*(sm) + (1-slope)*mean(bias) ;
% ma_sub = slope*(sm) + (1-slope)*bias ;

p.sub = ma_sub;
p.gm = ma_gm;

end
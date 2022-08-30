function [p, err] = determine_ma_pred_expt2a_V2(sm, bias, phi, slope)

if isempty(slope), slope = 0.5; end

%find the MA prediction based on individual subject data and the population average
ma_gm = ( slope*(sm-phi) + (1-slope)*(30+mean(bias)) );
ma_sub =( slope*(sm-phi) + (1-slope)*(30+bias) );

p.sub = ma_sub;
p.gm = ma_gm;

end
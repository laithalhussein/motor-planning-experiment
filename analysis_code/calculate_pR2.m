function [ pR2,p,s ] = calculate_pR2( original_data, model_data_small, model_data_large, p1, p2 )
%CALCULATE_pR2 ( original_data, model_data_small, model_data_large, p1, p2 )
% Put in data and the fits from both a small model and a large model in
% which the small model is embedded, along with the DOF of each (including
% offset).
% Output the pR^2 and p-value of the improvement.

if nargin<4
    computeSignificance = 0;
else
    computeSignificance = 1;
end

good_indices = ~(isnan(original_data)|isnan(model_data_small)|isnan(model_data_large));
original_data = original_data(good_indices);
model_data_small = model_data_small(good_indices);
model_data_large = model_data_large(good_indices);


rss1 = sum((original_data - model_data_small).^2);
rss2 = sum((original_data - model_data_large).^2);

pR2 = 1-rss2/rss1;

if(computeSignificance)
    n = length(original_data);
    F = (n-p2)/(p2-p1)*(1/(1-pR2)-1);
    p = 1-fcdf(F,p2-p1,n-p2);
else
    p = nan;
end

%return the F statistic and dof
s.F = F;
s.dof1 = p2-p1;
s.dof2 = n-p2;

end


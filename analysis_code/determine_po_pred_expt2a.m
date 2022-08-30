function p = determine_po_pred_expt2a(sm, phi, obs_ratio_all, obs_ratio_gm, slope)

if isempty(slope), slope = 0.5; end

%find the PO prediction based on individual subject data and the population average
obs_bound = phi;
tmp_phi = obs_bound *ones(size(sm));
new_theta = abs(sm)-tmp_phi;
new_theta(sm<=obs_bound) = 0;
po_sub = obs_ratio_all .* new_theta .* slope;
po_gm = obs_ratio_gm * new_theta * slope;

p.sub = po_sub;
p.gm = po_gm;

end
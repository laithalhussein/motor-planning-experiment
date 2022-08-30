function y = single_dof_model_flipped(b, x)
%in this function we set the offset of a linear regression model in
%accordance to the slope, i.e. we fit y = x*b + b, so there's only a single DOF

y = b*(x+-30) + (1-b)*30;

end
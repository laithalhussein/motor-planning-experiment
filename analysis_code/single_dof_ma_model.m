function y = single_dof_ma_model(b, x)
%in this function we set the offset of a linear regression model in
%accordance to the slope, i.e. we fit y = x*b + b, so there's only a single DOF
sm = x(:,1);
bias = x(:,2);
%y = b*(sm+-30) + (1-b)*(30+mean(bias));
%y = b*(sm) + (1-b)*mean(bias);
y = b*(sm) + (1-b).*(bias);
end
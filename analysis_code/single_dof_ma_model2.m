function y = single_dof_ma_model2(b, x)
%in this function we set the offset of a linear regression model in
%accordance to the slope, i.e. we fit y = x*b + b, so there's only a single DOF
sm = x(:,1);
bias = x(:,2);
y = b(1)*(sm+-30) + (1-b(1))*(30+mean(bias)) +b(2);

end
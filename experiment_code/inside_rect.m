function [obstacle_flag]=inside_rect(rect,xp,yp)
[A1,B1,C1]=get_edges(rect(2,:),rect(1,:));
[A2,B2,C2]=get_edges(rect(3,:),rect(2,:));
[A3,B3,C3]=get_edges(rect(4,:),rect(3,:));
[A4,B4,C4]=get_edges(rect(1,:),rect(4,:));

D1=A1*xp +B1*yp+C1;
D2=A2*xp +B2*yp+C2;
D3=A3*xp +B3*yp+C3;
D4=A4*xp +B4*yp+C4;

% if sign(D1)==sign(D2)==sign(D3)==sign(D4), obstacle_flag=1; else obstacle_flag=0; end
if numel(unique([sign(D1),sign(D2),sign(D3),sign(D4)]))==1, obstacle_flag=1; else obstacle_flag=0; end

return
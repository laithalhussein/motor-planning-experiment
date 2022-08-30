function [rect_points,circle_center,circle_center_flip]=make_obs_loc(cx,cy,theta,width,height)
%function to define the coordinates of a rotated rectangle, drawn to the
%screen as a polygon, NOT a texture
%this function will also return a 2nd point known to be at the center of
%one of the 2 short edges
theta=deg2rad(theta);
UL = [cx + (width/2) * cos(theta) - (height/2) * sin(theta) , cy + (height/2) * cos(theta) + (width/2) * sin(theta)];
UR = [cx - (width/2) * cos(theta) - (height/2) * sin(theta) , cy + (height/2) * cos(theta) - (width/2) * sin(theta)];
LL = [cx + (width/2) * cos(theta) + (height/2) * sin(theta) , cy - (height/2) * cos(theta) + (width/2) * sin(theta)];
LR = [cx - (width/2) * cos(theta) + (height/2) * sin(theta) , cy - (height/2) * cos(theta) - (width/2) * sin(theta)];

rect_points=[UL;UR;LR;LL]; %order is important here
%rect_points=[LL;LR;UR;UL];

circle_center=((UR-UL)/2)+UL; %check this
circle_center_flip=((LR-LL)/2) +LL;

end
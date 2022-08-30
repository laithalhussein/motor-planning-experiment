function [nx,ny]=make_circle(cx,cy,r)
%makes a circle centered at cx, cy with radius r
th = 0:pi/50:2*pi; %resolution -> ~0.001 degrees
nx= r * cos(th) + cx;
ny = r * sin(th) + cy;
return
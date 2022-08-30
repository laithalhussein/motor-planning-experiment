function [rx, ry] = rotatexy(x,y,phi)
% phi is in degrees
phi=phi*pi/180;
[theta r]=cart2pol(x,y);
[rx ry]=pol2cart(theta+phi,r);
return
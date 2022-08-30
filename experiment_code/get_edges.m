function [A,B,C]= get_edges(corner2,corner1)
%corners are passed as [x,y]
A=-(corner2(2)-corner1(2));
B=corner2(1)-corner1(1);
C=-(A*corner1(1)+B*corner1(2));

end
function [ V ] = xyz_1_To_XYZ_2( x,y,z,a,b,g,z0,y0,x0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
v0(:,1) = x;
v0(:,2) = y.*cos(g) - z.*sin(g);
v0(:,3) = y.*sin(g) + z.*cos(g);

v1(:,1) = v0(:,1).*cos(b) - v0(:,3).*sin(b);
v1(:,2) = v0(:,2);
v1(:,3) = v0(:,1).*sin(b) + v0(:,3).*cos(b);

V(:,1) = v1(:,1).*cos(a) - v1(:,2).*sin(a) + x0;
V(:,2) = v1(:,1).*sin(a) + v1(:,2).*cos(a) + y0;
V(:,3) = v1(:,3) + z0;

end


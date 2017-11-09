function [ V ] = XYZ_2_To_xyz_1( x,y,z,a,b,g,z0,y0,x0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
x1 = x - x0; 
y1 = y - y0;    
z1 = z - z0;

x2 = x1.*cos(a) + y1.*sin(a);
y2 = -x1.*sin(a) + y1.*cos(a); 
z2 = z1;

x3 = x2.*cos(b) + z2.*sin(b);
y3 = y2;  
z3 = -x2.*sin(b) + z2.*cos(b);

V(:,1) = x3;  
V(:,2) = y3.*cos(g) + z3.*sin(g);  
V(:,3) = -y3.*sin(g) + z3.*cos(g); 
end


function [ V ] = DPZ_To_XYZ( D,Ph, Z )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
try
    C = (D.^2-Z.^2).^(1/2);
    if (imag(C) ~= 0)
        throw(exception);
    end
catch 
    warning('Problem using function.  Assigning a value of 0.');
    C = 0;
end
V(:,1) = C.*cos(Ph);   
V(:,2) = -C.*sin(Ph);    
V(:,3) = Z;
end


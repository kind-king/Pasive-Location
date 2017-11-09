function [ W ] = xyz_To_DPZ( X,Y,Z )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
W(:,1) = (X.^2 + Y.^2 + Z.^2).^(1/2);
W(:,2) = ( atan(-Y./X) +  pi.*(X < 0) ) ./ pi .* 180;
W(:,3) = Z;


end


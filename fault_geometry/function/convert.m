function [anew] = convert(a,x0,y0)
%CONVERT Summary of this function goes here
%   Detailed explanation goes here
anew = a;
anew(:,1) = a(:,1) - x0;
anew(:,2) = a(:,2) - y0;
anew = anew / 1e3;
end


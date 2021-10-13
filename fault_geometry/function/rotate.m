function [b] = rotate(a,theta)
%ROTATE the coordiantes of fault geometry.
n = size(a,1);
m = size(a',1);
for i = 1: n
    b(i,2) = a(i,1)* sind(theta) + a(i,2) * cosd(theta);
    b(i,1) = a(i,1)* cosd(theta) - a(i,2) * sind(theta);
    if m == 3
    b(i,3) = a(i,3);
    end
end
end


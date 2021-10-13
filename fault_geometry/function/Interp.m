function [ycoor1]=a_Interp_Linear_Derivative(xcoor,ycoor,xcoor1,col,method)
% Jun 6, 2019

n=size(xcoor,1);
n1=size(xcoor1,1);
% Interpolation based on "method"
ycoor1 = interp1(xcoor, ycoor, xcoor1, method);


% figure(3)
% plot(xcoor,ycoor,'o','color',col);hold on;
% plot(xcoor1,ycoor1,'*','color',col);
% hold on;axis equal;

end
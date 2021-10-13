function [paleoloc,rates,loc_eqdyna] = observed_sliprates
%OBSERVED_SLIPRATES Summary of this function goes here
%   Detailed explanation goes here
% 06/16/2021.
% The data are from Dawson and Weldon, 2003.
% Geodetic locations are picked from the first data point.

addpath('E:\Backup\TAMU\Postdoc_work\Project#14-Cajon-Pass-Big-Bend-Terra\Smith_K_Loading\function');

x0 = 3.7e5; y0 = 3.8e6; % Origin
theta = 40;

paleoloc = [-160.542450889811,	-1.02422028868915; %BF
-88.5729384476460,	2.60096738540717;%MP
-68.4189296536326,	12.9274088985096;%FM
-42.1480708306533,	22.2660583874455;%3P
-20.4780566786978,	29.2188633408023;%EL
21.2535687819378,	39.6537262916120;%LR
28.7711882990793,	41.1736820727875;%PC
50.3325373551802,	46.6718474328023;%WW
71.6247010966469,	50.3032595021924;%LS
75.6720689211562,	50.6085392957337;%PT
105.622648060942,	56.1353768206584;%PL
98.8737529492582,	42.6180146755802;%Cn
125.113987371766,	40.5765816334258;%ML
172.007658926433,	39.2217647409170;]; %HL

% min, max, best, lon, lat;
rates = [31,37, 34, -119.8145, 35.2591; % Carrizo, 365, in TableB1
    %31, 37, 34, 100000, 100000; % Big Bend, 364
    20, 34, 27, -117.957, 34.487; % Mojave S, 374
    5, 20, 13, -116.896, 34.025; % San Bernardino S, 392
    2, 10, 6, -117.283, 34.059; % SJF, San Bernardino, 428 (Could be claremont?)
    6, 11, 8, -116.18, 33.282; % SJF, Clark. 419.
    ]; 
n = size(rates,1);

for i = 1:n
    if abs(rates(i,4))<361 
        [x,y,utmzone] = deg2utm(rates(i,5),rates(i,4)); 
        loc1(i,1) = x;
        loc1(i,2) = y;
    end
end

a = convert(loc1,x0,y0);
loc_eqdyna = rotate(a,theta);
end


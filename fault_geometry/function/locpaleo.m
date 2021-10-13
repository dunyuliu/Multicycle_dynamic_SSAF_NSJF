function [loc_eqdyna] = locpaleo(theta, x0, y0)
% The function is used to convert utm paleo sites coordinates to EQdyna xyx system.
% Documented on 07/26/2021.
% Last modified on 07/26/2021.

% The utm paleosites data come from Scharer and Yule (2020) updated table S1. 

loc=[         
3902410.1,   246358.9915;%BF
3858926.044, 303821.0623; %MP
3853881.81,	 325897.6374; %FM
3844149.048,	352025.051; %3P
3835545.989,	373094.4218; %EL
3816714.986,	411770.1022; %LR
3813047.107,	418505.9416;% PC
3803399.578,	438557.0458;%WW
3792495.062,	457202.0162;%LS
3790127.322, 460498.7099;%PT
3775109.264, 486994.7673;%PL
3762105.31213, 512926.570093; %Bu
3769092.47, 473136.0208; %Cn
3750661.744, 491924.9987; % ML
3719481.223, 526976.7757; %HL
];

loc1(:,1) = loc(:,2);
loc1(:,2) = loc(:,1);


a = convert(loc1,x0,y0);
loc_eqdyna = rotate(a,theta);

end

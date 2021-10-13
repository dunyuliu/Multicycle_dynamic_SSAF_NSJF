clear all; close all;
% Documented on 07/26/2021.
% Last updated on 09/12/2021.
% This script is to convert CFM5.2 with UTM coordinates to EQdyna xyz.
% 
% in /function/locpaleo.m add Bu site UTM locations.
% Paleo_sites_loc_inEQdyna.txt is saved from ploc array. 
% To visualize and preprocess the Cajon Pass fault trace for meshing.
addpath('./function/')
addpath('./rawfaulttracedata/') % stores the raw CFM5.2 data with UTM coords.
a10 = load('PNRA-SJFZ-SBRN-SJ-Claremont-coor.txt');
a11 = load('PNRA-SJFZ-SJVC-SJ-Claremont-coor.txt');
a12 = load('PNRA-SJFZ-ANNA-Clark-combined.txt');
a13 = load('PNRA-SJFZ-SJVC-Claremont-coor.txt');
a20 = load('SAFS-SAFZ-CLCZ-CHLM-coor.txt');
a30 = load('SAFS-SAFZ-CLCZ-CRRZ-coor.txt');
a40 = load('SAFS-SAFZ-MJVS-coor.txt');
a50 = load('SAFS-SAFZ-MULT1-coor.txt');
a60 = load('SAFS-SAFZ-MULT2-coor.txt');
a70 = load('SAFS-SAFZ-PARK-coor.txt');
a80 = load('SAFS-SAFZ-SBMT-MillC-coor.txt');
a90 = load('SAFS-SAFZ-SBMT-MissionC-coor.txt');
a100 = load('SAFS-SAFZ-SBMT-SAF-coor.txt'); 
%-------------------Parameters--------------------------------------------%
x0 = 3.7e5; y0 = 3.8e6; % Origin
x1 = [-160, 149.9]; x2 = [-113.8, 91.38]; % Two points handpicked on the left 
% part of the big bend. The angle, 51.some, is used to rotate the
% cooridantes.
d0 = -3.5; d1 = -2.5; % d0 and d1, in km, are lower and upper cap of zcoors to retrieve fault traces.
% d0 = -10.5; d1 = -9.5; % d0 and d1, in km, are lower and upper cap of zcoors to retrieve fault traces.
dx = 0.4; %This is actually a spatial filter, it does not necessarily coincide with 
%the dx in meshing.
xx =  -350: dx : 250; nx = size(xx',1);
method = 'Linear'

%-------------------------------------------------------------------------%
a1 = convert(a10,x0,y0); %SB
b1 = convert(a11,x0,y0);
b2 = convert(a12,x0,y0);
b3 = convert(a13,x0,y0);
a2 = convert(a20,x0,y0); %CHLM
a3 = convert(a30,x0,y0); %CRRZ
a4 = convert(a40,x0,y0); %MJVS
a5 = convert(a50,x0,y0);
a6 = convert(a60,x0,y0);
a7 = convert(a70,x0,y0);
a8 = convert(a80,x0,y0);
a9 = convert(a90,x0,y0);
a10 = convert(a100,x0,y0); %SBMT
 
theta = atand((x2(2) - x1(2))/(x2(1) - x1(1)))
theta = abs(theta)
theta = 40

[ploc] = locpaleo(theta, x0, y0);

asjf = [a1;b1];
asjfclare = b3;
asjfanna = b2;
asaf = [a2;a3;a4;a10];%
% numbers of controlling points in each segment are [402;1199;664;688;]
% 

asjf1 = rotate(asjf,theta);
asjf2 = rotate(asjfclare,theta);
asjf3 = rotate(asjfanna,theta);
asaf1 = rotate(asaf,theta);

figure(1)
plot3(a1(:,1),a1(:,2),a1(:,3),'r*'); hold on;
plot3(b1(:,1),b1(:,2),b1(:,3),'ko'); hold on;
plot3(b2(:,1),b2(:,2),b2(:,3),'b*'); hold on;
plot3(b3(:,1),b3(:,2),b3(:,3),'r*'); hold on;
plot3(a2(:,1),a2(:,2),a2(:,3),'b*'); hold on;
plot3(a3(:,1),a3(:,2),a3(:,3),'k*'); hold on;
plot3(a4(:,1),a4(:,2),a4(:,3),'m*'); hold on;
plot3(a5(:,1),a5(:,2),a5(:,3),'c*'); hold on;
plot3(a6(:,1),a6(:,2),a6(:,3),'ko'); hold on;
plot3(a7(:,1),a7(:,2),a7(:,3),'bo'); hold on;
plot3(a8(:,1),a8(:,2),a8(:,3),'ko'); hold on;
plot3(a9(:,1),a9(:,2),a9(:,3),'g*'); hold on;
plot3(a10(:,1),a10(:,2),a10(:,3),'y*'); hold on;
set(gcf,'color','white');
figure(2)
plot3(a1(:,1),a1(:,2),a1(:,3),'r*'); hold on;
plot3(b1(:,1),b1(:,2),b1(:,3),'ko'); hold on;
plot3(b2(:,1),b2(:,2),b2(:,3),'b*'); hold on;
plot3(b3(:,1),b3(:,2),b3(:,3),'r*'); hold on;
plot3(a2(:,1),a2(:,2),a2(:,3),'b*'); hold on;
plot3(a3(:,1),a3(:,2),a3(:,3),'k*'); hold on;
plot3(a4(:,1),a4(:,2),a4(:,3),'m*'); hold on;
plot3(a5(:,1),a5(:,2),a5(:,3),'c*'); hold on;
%plot3(a6(:,1),a6(:,2),a6(:,3),'ko'); hold on;
plot3(a7(:,1),a7(:,2),a7(:,3),'bo'); hold on;
%plot3(a8(:,1),a8(:,2),a8(:,3),'ko'); hold on;
%plot3(a9(:,1),a9(:,2),a9(:,3),'g*'); hold on;
plot3(a10(:,1),a10(:,2),a10(:,3),'g*'); hold on;

t1 = findtrace(rotate(a1,theta),d0,d1); ends(1,1) = min(t1(:,1));ends(1,2) = max(t1(:,1));
p1 = findtrace(rotate(b1,theta),d0,d1); ends(1,1) = min(p1(:,1));ends(1,2) = max(p1(:,1));
p2 = findtrace(rotate(b2,theta),d0,d1); ends(1,1) = min(p2(:,1));ends(1,2) = max(p2(:,1));
p3 = findtrace(rotate(b2,theta),d0,d1); ends(1,1) = min(p3(:,1));ends(1,2) = max(p3(:,1));
t2 = findtrace(rotate(a2,theta),d0,d1); ends(2,1) = min(t2(:,1));ends(2,2) = max(t2(:,1));
t3 = findtrace(rotate(a3,theta),d0,d1); ends(3,1) = min(t3(:,1));ends(3,2) = max(t3(:,1));
t4 = findtrace(rotate(a4,theta),d0,d1); ends(4,1) = min(t4(:,1));ends(4,2) = max(t4(:,1));
t5 = findtrace(rotate(a5,theta),d0,d1); ends(5,1) = min(t5(:,1));ends(5,2) = max(t5(:,1));
%t6 = findtrace(a6);
t7 = findtrace(rotate(a7,theta),d0,d1); ends(6,1) = min(t7(:,1));ends(6,2) = max(t7(:,1));
%t8 = findtrace(a8);
%t9 = findtrace(a9);
t10 = findtrace(rotate(a10,theta),d0,d1); ends(7,1) = min(t10(:,1));ends(7,2) = max(t10(:,1));

tsjf0 = findtrace(asjf1,d0,d1);
tsjf1 = findtrace(asjf2,d0,d1);
tsjf2 = findtrace(asjf3,d0,d1);
tsaf0 = findtrace(asaf1,d0,d1);
% tsjf and tsaf are controlling points in Cartesian coordiantes 
% before intepolation. tsjf is in flt1coor.txt. tsaf is in flt2coor.txt.
tsjf = unique(tsjf0(:,1:2),'rows');
tsjf1 = unique(tsjf1(:,1:2),'rows');
tsjf2 = unique(tsjf2(:,1:2),'rows');
tsaf = unique(tsaf0(:,1:2),'rows');
% %axis equal;
% %%SJF
% nf = size(tsjf,1);
% nxtmp1=0;
% for i=1:nx
%     for j=1:nf
%         if xx(i)> (tsjf(1,1)) && xx(i) < (tsjf(nf,1))
%             nxtmp1 = nxtmp1 + 1;
%             xtmp1(nxtmp1,1) = xx(i);
%             break;
%         end
%     end    
% end
% [ycoor1] = Interp(tsjf(:,1),tsjf(:,2),xtmp1,'k',method);
% %%SAF
% nf = size(tsaf,1);
% nxtmp2=0;
% for i=1:nx
%     for j=1:nf
%         if xx(i)> (tsaf(1,1)) && xx(i) < (tsaf(nf,1))
%             nxtmp2 = nxtmp2 + 1;
%             xtmp2(nxtmp2,1) = xx(i);
%             break;
%         end
%     end    
% end
% [ycoor2] = Interp(tsaf(:,1),tsaf(:,2),xtmp2,'k',method);
% xy1(:,1)=xtmp1;xy1(:,2)=ycoor1;
% xy2(:,1)=xtmp2;xy2(:,2)=ycoor2;

%Plotting
figure(3)
plot(t1(:, 1), t1(:,2), 'r*'); hold on;
plot(t2(:, 1), t2(:,2), 'b*'); hold on;
plot(t3(:, 1), t3(:,2), 'k*'); hold on;
plot(t4(:, 1), t4(:,2), 'm*'); hold on;
plot(t5(:, 1), t5(:,2), 'c*'); hold on;
%plot(t6(:, 1), t6(:,2), 'ro'); hold on;
plot(t7(:, 1), t7(:,2), 'bo'); hold on;
%plot(t8(:, 1), t8(:,2), 'ko'); hold on;
%plot(t9(:, 1), t9(:,2), 'mo'); hold on;
plot(t10(:, 1), t10(:,2), 'g*'); hold on;
axis equal;


h = figure(4);
set(h, 'Position', [100 100 1500 450]);
plot(tsjf(:, 1), tsjf(:,2), 'r*'); hold on;
plot(tsjf1(:, 1), tsjf1(:,2), 'k*'); hold on;
%plot(tsjf2(:, 1), tsjf2(:,2), 'g*'); hold on;
plot(tsaf(:, 1), tsaf(:,2), 'b*'); hold on;
plot(ploc(:,1),ploc(:,2),'k*','markersize',25);

axis equal;
xlabel('NW-SE (km)'); ylabel('SW-NE (km)'); 
ylim([-30 100]);xlim([-250 170]);
set(gcf,'color','white');
set(gca, 'fontsize',14, 'fontweight','bold');

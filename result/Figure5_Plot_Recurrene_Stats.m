clear all; close all;
% Created on 07/27/2021.
% Last updated on 10/13/2021.
% Author: Dunyu Liu, dliu@ig.utexas.edu.

% Parameters:
% -slipthreshod: at a site, over which value the earthquake is counted.
% -mod: model id. 1-3 for Model B, A, C, respectively.
% -nsite: which paleo site to pick. Currently available options are
% 1,3,5,8.
%   --1: BF;
%   --3: FM;
%   --5: EL;
%   --8: WW:
% [Note: the nsmp id of paleo sites are hand picked from
% Figure1_Prep_SAF_SJF_Trace.m. The UTM locations of various sites are in
% locpaleo.m in ./fault_geometry/function/.

% Key arrays:
% -plocall is n by 2 array. The first column stores the integer index of a
% site in the EQdyna nsmp master-slave node file. 
% The second column ranges from 1-3, indicating the fault id. 1 for SAF. 2
% for Claremont. 3 for Clark.

% Dependency:
% -Paleo_sites_loc_inEQdyna.txt stores the EQdyna_xy locations of the paleo sites.
%  There notations

% With the current parameter set, BF from Model A will be plotted.
mod = 2;
slipthreshold = 0.5;
nsite = 1;

plocall = [358,3;%BF 1 
    33,3;%MP 2 
    818,3;%FM 3
    44,3;%3P 4 
    1058,3;%EL 5
    57,3;%LR 6
    59,3;%PC 7
    1412,3;%WW 8
    69,3;%LS 9
    760,3; % 12 km north of FM; 10
    706,3; % 22 km north of FM; 11
    ];
ploc(1:2) = plocall(nsite,1:2);

lockdepth = 14e3; rig = 3464*3464*2670; alpha = 1;
if mod == 1
    path = './work_vis4.2_fs0.3/'; % Model B
    filenum = 4;
    filetag = [1,154,1378,2888;];
    icstart = 1;
elseif mod == 2
    path = './work_vis7_fs0.5/'; % Model A
    filenum = 3;
    filetag = [1,1256,2929;];
    icstart = 1;
elseif mod == 3
    path = './work_vis12_fs0.7/'; % Model C
    filenum = 3;
    filetag = [1,1436,3284;]; 
    icstart = 1;
end
if isfile(strcat(path,'data.mat'))
    data = load(strcat(path,'data.mat'));
    ic = data.ic; res = data.res; tinte = data.tinte;
else
    k = 1;
    for kk = 1: filenum
        if kk == 1
            ic = load(strcat(path,'cyclelog.txt',num2str(filetag(k,kk))'));
            res = load(strcat(path,'totalop.txt',num2str(filetag(k,kk))'));
            tinte  = load(strcat(path,'interval.txt',num2str(filetag(k,kk))'));
        else
            k
            ictmp = load(strcat(path,'cyclelog.txt',num2str(filetag(k,kk))));
            restmp = load(strcat(path,'totalop.txt',num2str(filetag(k,kk))));
            tintetmp  = load(strcat(path,'interval.txt',num2str(filetag(k,kk))));
            ic(1,2) = ictmp(1,2);
            res = [res;restmp;];
            tinte = [tinte;tintetmp;];
        end
    end
    %saveas(strcat(path,'data.mat'),'ic','res','tinte');
end
nt = size(tinte,1);
ttot(1,1)=0;
for i = 1:nt-1
    ttot(i+1,1)=ttot(i,1)+tinte(i+1,1)/1e3;
end

if exist('mesh.mat','file')
    C = load('mesh.mat');
    vert2 = C.vert2;
    nsmpnv2 = C.nsmpnv2;
    nsmp2 = C.nsmp2;
else
    vert2 = load(strcat('mesh/vert.txt')); vert2 = vert2/1e3;
    nsmpnv2 = load(strcat('mesh/nsmpnv.txt'));
    nsmp2 = load(strcat('mesh/nsmp.txt'));
end

totft = 3;
nft = [295,178,1769]; maxftnode = 1769;
tag = nft(1); nftsum(1) = tag;
for i = 2:totft
    tag = tag + nft(i);
    nftsum(i) = tag; 
end
ntotnd = sum(nft)
x1 = vert2(nsmp2(1:nft(1),1),:);
x2 = vert2(nsmp2(maxftnode+1: maxftnode + nft(2),1),:);
x3 = vert2(nsmp2(maxftnode*2+1: maxftnode*2+ nft(3),1),:);
nx1 = size(x1,1)
nx2 = size(x2,1)
nx3 = size(x3,1)
col = ['b','r','c','m','k','g','b'];
col2 = ['b:','r:','c:','m:','k','g','b'];


np = size(ploc,1);
for i = 1: np
    if ploc(i,2)>1
        plocnew(i) = nftsum(ploc(i,2)-1) + ploc(i,1);
    else
        plocnew(i) = ploc(i,1);
    end
end

ntag = 0;
ictag = 0
mom = 0;
ic(1) = icstart;
for i = ic(1): ic(2)
    ictag = ictag + 1;
    i;
    tmp = res((ictag-1)*ntotnd+1:ictag*ntotnd, :);
    ns(i,:) = tmp(1:ntotnd, 2);
    ss(i,:) = tmp(1:ntotnd, 1);
    slip(i,:) = tmp(1:ntotnd, 3);
    rupt(i,:) = tmp(1:ntotnd, 5);
    maxslip = max(slip(i,765:1728));
    for jj = 1:nftsum(1)
        if jj ==1 
            mom = mom + slip(i,jj)* 0.5 * ((x1(1,1)-x1(2,1))^2 + (x1(1,2)-x1(2,2))^2)^0.5 * rig;
        elseif jj>1&& jj<nftsum(1)
            mom = mom + slip(i,jj)* 0.5* rig * (((x1(jj-1,1)-x1(jj,1))^2 + (x1(jj-1,2)-x1(jj,2))^2)^0.5 + ((x1(jj,1)-x1(jj+1,1))^2 + (x1(jj,2)-x1(jj+1,2))^2)^0.5);
        elseif jj == nftsum(1)
            mom = mom + slip(i,jj)*0.5*rig* ((x1(jj-1,1)-x1(jj,1))^2 + (x1(jj-1,2)-x1(jj,2))^2)^0.5;
        end
    end
    for jj = 1+nftsum(1):nftsum(2)
        if jj ==1+nftsum(1) 
            mom = mom + slip(i,jj)* 0.5 * ((x2(jj-nftsum(1),1)-x2(jj-nftsum(1)+1,1))^2 + (x2(jj-nftsum(1),2)-x2(jj-nftsum(1)+1,2))^2)^0.5 * rig;
        elseif jj>1+nftsum(1)&& jj<nftsum(2)
            mom = mom + slip(i,jj)* 0.5* rig * (((x2(jj-nftsum(1)-1,1)-x2(jj-nftsum(1),1))^2 + (x2(jj-nftsum(1)-1,2)-x2(jj-nftsum(1),2))^2)^0.5 + ((x2(jj-nftsum(1),1)-x2(jj-nftsum(1)+1,1))^2 + (x2(jj-nftsum(1),2)-x2(jj-nftsum(1)+1,2))^2)^0.5);
        elseif jj == nftsum(2)
            mom = mom + slip(i,jj)*0.5*rig* ((x2(jj-nftsum(1)-1,1)-x2(jj-nftsum(1),1))^2 + (x2(jj-nftsum(1)-1,2)-x2(jj-nftsum(1),2))^2)^0.5;
        end
    end
    mom = alpha * mom * lockdepth*1e3;
    mag = 2/3*log10(mom)-6.07;    

    if slip(i,plocnew)>slipthreshold
            i
        ntag = ntag + 1;
        icyc(ntag) = i;
        timestamp(ntag) = ttot(i);
        magrec(ntag) = mag;
        maxsliprec(ntag) = maxslip;
        if ntag > 1
            recurrence(ntag-1) = timestamp(ntag) - timestamp(ntag-1);
        end
    end
end
meanrecurr = mean(recurrence)
stand = std(recurrence)
COV = stand/meanrecurr

meanslip = mean(maxsliprec)
stdslip = std(maxsliprec)

recur_bin = 0:0.05:0.6; nm = size(recur_bin',1);
recur_binplot = 0.025:0.05:0.575;
count(1:nm-1) = 0;

for i = 1: length(recurrence)
    for k = 2:nm
        if (recurrence(i)>=recur_bin(k-1) && recurrence(i)<recur_bin(k))
            count(k-1) = count(k-1) + 1;
        end
    end
end

position=[100 100 750 350];
set(0,'DefaultFigurePosition', position);
fig1=figure(1);
h1 = subplot(1,2,1);
bar(recurrence); hold on;

xlabel('The Number of Interval from Past to Present');
ylabel('Recurrence Interval (kyrs)');
set(gcf, 'color', 'white');
text(h1,0.025,0.9,'(e)','Units','normalized','FontSize',12);



h2 = subplot(1,2,2);
sumcount = sum(count);
yyaxis left;
bar(recur_binplot,count); hold on;
ylabel('Number of Intervals');

text(0.025,0.9,'(f)','Units','normalized','FontSize',12);

yyaxis right;
cdfplot(recurrence);hold on;
ylabel('Frequency');xlabel('Recurrence Intervals (kyrs)');
title('');
text(h2,0.025,0.6, strcat('COV=',num2str(COV,'%4.2f')),'Units','normalized','FontSize',12);
text(h2,0.025,0.55, strcat('Mean=',num2str(meanrecurr,'%4.2f')),'Units','normalized','FontSize',12);
text(h2,0.025,0.5, strcat('Std=',num2str(stand,'%4.2f')),'Units','normalized','FontSize',12);
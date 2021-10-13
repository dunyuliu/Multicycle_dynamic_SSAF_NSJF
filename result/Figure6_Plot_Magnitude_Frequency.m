clear all; close all;
% Created on 9/13/2021.
% Last modified on 10/13/2021.
% Author: Dunyu Liu, dliu@ig.utexas.edu.
% The script is to plot Figure 6 in the JGR manuscript.
slipthreshold = 0.5;
m = 2;

lockdepth = 22e3; rig = 3500*3500*3000; alpha = 1;
if m == 1
%     path = './work_vis4.2_fs0.3/'
%     filenum = 4;
%     filetag = [1,154,1378,2888;];
elseif m == 2
    path = './work_vis7_fs0.5/'
    filenum = 3;
    filetag = [1,1256,2929;];
end

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

ntag = 0;
ictag = 0

for i = ic(1): ic(2)
    mom = 0;
    len = 0;
    ictag = ictag + 1
    i
    tmp = res((ictag-1)*ntotnd+1:ictag*ntotnd, :);
    ns(i,:) = tmp(1:ntotnd, 2);
    ss(i,:) = tmp(1:ntotnd, 1);
    slip(i,:) = tmp(1:ntotnd, 3);
    rupt(i,:) = tmp(1:ntotnd, 5);
    maxslip = max(slip(i,:));
    for jj = 1:nftsum(1)
        if jj ==1 
            mom = mom + slip(i,jj)* 0.5 * ((x1(1,1)-x1(2,1))^2 + (x1(1,2)-x1(2,2))^2)^0.5 * rig;
        elseif jj>1&& jj<nftsum(1)
            mom = mom + slip(i,jj)* 0.5* rig * (((x1(jj-1,1)-x1(jj,1))^2 + (x1(jj-1,2)-x1(jj,2))^2)^0.5 + ((x1(jj,1)-x1(jj+1,1))^2 + (x1(jj,2)-x1(jj+1,2))^2)^0.5);
        elseif jj == nftsum(1)
            mom = mom + slip(i,jj)*0.5*rig* ((x1(jj-1,1)-x1(jj,1))^2 + (x1(jj-1,2)-x1(jj,2))^2)^0.5;
        end
        if (slip(i,jj)>0) 
            if jj ==1 
                len = len + 0.5 * ((x1(1,1)-x1(2,1))^2 + (x1(1,2)-x1(2,2))^2)^0.5; 
            elseif jj>1&& jj<nftsum(1)
                len = len + 0.5* (((x1(jj-1,1)-x1(jj,1))^2 + (x1(jj-1,2)-x1(jj,2))^2)^0.5 + ((x1(jj,1)-x1(jj+1,1))^2 + (x1(jj,2)-x1(jj+1,2))^2)^0.5);
            elseif jj == nftsum(1)
                len = len + 0.5* ((x1(jj-1,1)-x1(jj,1))^2 + (x1(jj-1,2)-x1(jj,2))^2)^0.5;
            end
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
        if slip(i,jj) > 0
            if jj ==1+nftsum(1) 
                len = len + 0.5 * ((x2(jj-nftsum(1),1)-x2(jj-nftsum(1)+1,1))^2 + (x2(jj-nftsum(1),2)-x2(jj-nftsum(1)+1,2))^2)^0.5;
            elseif jj>1+nftsum(1)&& jj<nftsum(2)
                len = len + 0.5* (((x2(jj-nftsum(1)-1,1)-x2(jj-nftsum(1),1))^2 + (x2(jj-nftsum(1)-1,2)-x2(jj-nftsum(1),2))^2)^0.5 + ((x2(jj-nftsum(1),1)-x2(jj-nftsum(1)+1,1))^2 + (x2(jj-nftsum(1),2)-x2(jj-nftsum(1)+1,2))^2)^0.5);
            elseif jj == nftsum(2)
                len = len + 0.5* ((x2(jj-nftsum(1)-1,1)-x2(jj-nftsum(1),1))^2 + (x2(jj-nftsum(1)-1,2)-x2(jj-nftsum(1),2))^2)^0.5;
            end                
        end
    end
    for jj = 1+nftsum(2):nftsum(3)
        if jj ==1+nftsum(2) 
            mom = mom + slip(i,jj)* 0.5 * ((x3(jj-nftsum(2),1)-x3(jj-nftsum(2)+1,1))^2 + (x3(jj-nftsum(2),2)-x3(jj-nftsum(2)+1,2))^2)^0.5 * rig;
        elseif jj>1+nftsum(2)&& jj<nftsum(3)
            mom = mom + slip(i,jj)* 0.5* rig * (((x3(jj-nftsum(2)-1,1)-x3(jj-nftsum(2),1))^2 + (x3(jj-nftsum(2)-1,2)-x3(jj-nftsum(2),2))^2)^0.5 + ((x3(jj-nftsum(2),1)-x3(jj-nftsum(2)+1,1))^2 + (x3(jj-nftsum(2),2)-x3(jj-nftsum(2)+1,2))^2)^0.5);
        elseif jj == nftsum(3)
            mom = mom + slip(i,jj)*0.5*rig* ((x3(jj-nftsum(2)-1,1)-x3(jj-nftsum(2),1))^2 + (x3(jj-nftsum(2)-1,2)-x3(jj-nftsum(2),2))^2)^0.5;
        end
        if slip(i,jj) > 0
            if jj ==1+nftsum(2) 
                len = len + 0.5* ((x3(jj-nftsum(2),1)-x3(jj-nftsum(2)+1,1))^2 + (x3(jj-nftsum(2),2)-x3(jj-nftsum(2)+1,2))^2)^0.5;
            elseif jj>1+nftsum(2)&& jj<nftsum(3)
                len = len + 0.5* (((x3(jj-nftsum(2)-1,1)-x3(jj-nftsum(2),1))^2 + (x3(jj-nftsum(2)-1,2)-x3(jj-nftsum(2),2))^2)^0.5 + ((x3(jj-nftsum(2),1)-x3(jj-nftsum(2)+1,1))^2 + (x3(jj-nftsum(2),2)-x3(jj-nftsum(2)+1,2))^2)^0.5);
            elseif jj == nftsum(3)
                len = len + 0.5* ((x3(jj-nftsum(2)-1,1)-x3(jj-nftsum(2),1))^2 + (x3(jj-nftsum(2)-1,2)-x3(jj-nftsum(2),2))^2)^0.5;
            end        
        end
    end    
    mom = alpha * mom * min(len*1e3,lockdepth)*1e3;
    mag = 2/3*log10(mom)-6.07;    

    icyc(i) = i;
    timestamp(i) = ttot(i);
    magrec(i) = mag;
    mo(i) = mom;
    if i == 1
        moaccum(i) = mom;
    else
        moaccum(i) = moaccum(i-1) + mom;
    end
    maxsliprec(i) = maxslip;
    if i > 1
        recurrence(i-1) = timestamp(i) - timestamp(i-1);
    end
end
meanrecurr = mean(recurrence)
stand = std(recurrence)
COV = stand/meanrecurr

meanslip = mean(maxsliprec)
stdslip = std(maxsliprec)

position=[100 100 550 350];
set(0,'DefaultFigurePosition', position);
fig1=figure(1);
bar(recurrence);
text(floor(ntag/2),max(recurrence)*1.3, strcat('COV=',num2str(COV)));
text(floor(ntag/2),max(recurrence)*1.25, strcat('Mean=',num2str(meanrecurr)));
text(floor(ntag/2),max(recurrence)*1.2, strcat('Std=',num2str(stand)));
ylim([0 max(recurrence)*1.3])
xlabel('The Number of Interval from Past to Present');
ylabel('Recurrence Interval (kyrs)');
set(gcf, 'color', 'white');
%set(gca, 'fontsize', 12, 'fontweight', 'bold');

fig1=figure(2);
magrange = 5.5:0.1:8.5; nm = size(magrange',1);
count(1:nm-1) = 0;
magrangeplot = 5.55:0.1:8.45;
for i = 1: ic(2)
    for k = 2:nm
        if (magrec(i)>=magrange(k-1) && magrec(i)<magrange(k))
            count(k-1) = count(k-1) + 1;
        end
    end
end
bar(magrangeplot,count);
xlabel('Magnitudes');
ylabel('The Number of Earthquakes');

crat = count/ic(2);
crat = log10(crat);

figure(3)
plot(magrangeplot, crat,'r*'); hold on;
plot(magrangeplot, -magrangeplot,'bo');
legend('log10N/Ntot','-Mag');
xlabel('-Magnitudes');
ylabel('log10N/Ntot or - Mag');

XX = [ones(length(-magrangeplot),1) -magrangeplot'];
b = XX\crat'

h = figure(4);
set(h,'Position',[50 50 1000 400]);
subplot(1,3,1)
plot(timestamp, moaccum);
xlabel('Thousand Years');
ylabel('Cummulative Moment Release (Nm)')
set(gca,'Fontsize',10,'Fontweight','bold');
text(0.025,0.9,'(a)','Units','normalized','FontSize',12);
set(gca,'Fontsize',10,'Fontweight','bold');

sumcount = sum(count);
subplot(1,3,2)
yyaxis left;
bar(magrangeplot,count); 
ylabel('Number of Events');

yyaxis right;
cdfplot(magrec);hold on;
xlim([5 8.5]);
title('');
xlabel('Magnitude');
ylabel('Freqeuncy');
text(0.025,0.9,'(b)','Units','normalized','FontSize',12);
set(gca,'Fontsize',10,'Fontweight','bold');

ntag = 0;
for i = 1:ic(2)
    if magrec(i)>6.6
        ntag = ntag + 1;
        magrec1(ntag) = magrec(i);
    end
end

magrange1 = 6.5:0.1:8.5; nm1 = size(magrange1',1);
magrangeplot1 = 6.55:0.1:8.45 ; magrangeplot1 = magrangeplot1-0.03;
count1(1:nm1-1) = 0;

for i = 1: ntag
    for k = 2:nm1
        if (magrec1(i)>=magrange1(k-1) && magrec1(i)<magrange1(k))
            count1(k-1) = count1(k-1) + 1;
        end
    end
end

sumcount1 = sum(count1);


mag_obs = 6.51:0.1:8.01;
mag_obs_plot1 = 6.52:0.1:8.02;
mag_obs_plot2 = mag_obs_plot1 + 0.01;

count_obs1 = [1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,0]; 
count_obs2 = [0,0,0,0,0,2,7,3,4,7,4,2,4,0,1,0];
magrec_obs1 = [6.4, 6.5, 6.9, 7.2, 7.3, 7.5, 7.8]; % Observation for historic events. From Scharer and Yule (2020).
% obs2: Observation for paleo events. From Scharer and Yule (2020).
% 10/13/2021. Recontructed following Kate's 10/09/2021 email figure.
magrec_obs2 = [7.0, 7.0, ...
    7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1, ...
    7.2, 7.2, 7.2, ...
    7.3, 7.3, 7.3, 7.3, ...
    7.4, 7.4, 7.4, 7.4, 7.4, 7.4, 7.4, ...
    7.5, 7.5, 7.5, 7.5, ...
    7.6, 7.6, ...
    7.7, 7.7, 7.7, 7.7, ...
    7.9];
subplot(1,3,3)
yyaxis left; 
h5=bar(magrangeplot1,count1, 'b');hold on;h5.FaceAlpha=0.3;
h6=bar(mag_obs_plot1,count_obs2, 'r');hold on;h6.FaceAlpha=0.3;
h7=bar(mag_obs_plot2,count_obs1, 'k');hold on;h7.FaceAlpha=0.3;
ylabel('Number of events');

yyaxis right;
h1=cdfplot(magrec1);hold on;  
h2=cdfplot(magrec_obs2);hold on;
h3=cdfplot(magrec_obs1);hold on;
set(h1,'color','b' ,'Linewidth', 2);
set(h2,'color','r','Linewidth', 2);
set(h3,'color','k','Linewidth', 2);
xlim([6.0 8.5]);
title('');
xlabel('Magnitude');
ylabel('Freqeuncy');
text(0.025,0.9,'(c)','Units','normalized','FontSize',12);
set(gca,'Fontsize',10,'Fontweight','bold');
set(gcf, 'color', 'white');


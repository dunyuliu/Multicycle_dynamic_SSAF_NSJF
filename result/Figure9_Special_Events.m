clear all; close all;
% Last modified on 10/13/2021.
% The script is used to plot Figure 9 in the JGR manuscript. 
% The figure shows representitive and special events in Model A's sequence.
% Dependency:
% -mesh.mat or ./mesh/
% Locking depth is assigned to be 22 km to reproduce the magnitude 7.9 of
% the 1857 Fort Tejon earthquake.

nfig = 10;

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
ntotnd = sum(nft);
x1 = vert2(nsmp2(1:nft(1),1),:);
x2 = vert2(nsmp2(maxftnode+1: maxftnode + nft(2),1),:);
x3 = vert2(nsmp2(maxftnode*2+1: maxftnode*2+ nft(3),1),:);

%for k = 1:3
k = 2;
model = k;
lockdepth = 22e3; rig = 3500*3500*3000; alp = 1;
clear accumulateslip ttot ic tinte res;
%%
if model == 1 
    path = './work_vis4.2_fs0.3/'
    mius = 0.3;
    cycstart = 1;%The first rupture of most of the system.
    Mname = 'Model B'%'42fs03'
    k1 = 1;
    filetag = [1,154,1378,2888;];
    filenum = 4;       
elseif model == 2
    path = './work_vis7_fs0.5/'
    mius = 0.5;
    cycstart = 1;
    Mname = 'Model A'%'7fs05'
     k1 = 1;
    filetag = [1,1256,2929;];
    filenum = 3;         
elseif model == 3
    path = './work_vis12_fs0.7/'
    mius = 0.7;
    cycstart = 1;  
    Mname = 'Model C'%'12fs07'
      k1 = 1;
    filetag = [1,1436,3284;];
    filenum = 3;         
end

nx1 = size(x1,1)
nx2 = size(x2,1)
nx3 = size(x3,1)
col = ['b','r','k','m','g','c','y'];
col2 = ['b:','r:','k:','m:','g','c','y'];

for kk = 1: filenum
    if kk == 1
        ic = load(strcat(path,'cyclelog.txt',num2str(filetag(k1,kk))));
        res = load(strcat(path,'totalop.txt',num2str(filetag(k1,kk))));
        tinte  = load(strcat(path,'interval.txt',num2str(filetag(k1,kk))));
    else
        ictmp = load(strcat(path,'cyclelog.txt',num2str(filetag(k1,kk))));
        restmp = load(strcat(path,'totalop.txt',num2str(filetag(k1,kk))));
        tintetmp  = load(strcat(path,'interval.txt',num2str(filetag(k1,kk))));
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

ntag = 0;
ictag = 0
accumulateslip = zeros(ic(2) - ic(1) + 1,ntotnd);
tspan = ttot(ic(2)-ic(1)+1) - ttot(1);    
for i = ic(1): ic(2)
    ictag = ictag + 1;
    tmp = res((ictag-1)*ntotnd+1:ictag*ntotnd, :);
    ns(i,:) = tmp(1:ntotnd, 2);
    ss(i,:) = tmp(1:ntotnd, 1);
    slip(i,:) = tmp(1:ntotnd, 3);
    rupt(i,:) = tmp(1:ntotnd, 5);
    for j = 1:length(rupt(i,:))
        if rupt(i,j)>500
            rupt(i,j) = NaN;
        end
    end
    maxslip = max(slip(i,:));

    if i == 856 || i ==1041 || i == 2047 || i == 355 || i==2704 || i==3677 || i==3904 || i==3909
        i
        ntag = ntag + 1;
        h = figure(1);
        set(h, 'Position', [50 50 800 1000]);
        if i == 856
            subplot(nfig,1,1)
            mtag = 1;
        elseif i ==1041
            subplot(nfig,1,2)
            mtag = 2;
        elseif i == 2047
            subplot(nfig,1,3)
            mtag = 3;
        elseif i == 355
            subplot(nfig,1,4)
            mtag = 4;
        elseif i == 2704
            subplot(nfig,1,5)
            mtag = 5;
        elseif i == 3677 
            subplot(nfig,1,6)
            mtag = 6;
        elseif i == 3904
            subplot(nfig,1,7)
            mtag = 7;
        elseif i == 3909
            subplot(nfig,1,8)
            mtag = 8;
        end
        yyaxis left;
        plot(x1(:,1),slip(i,1:nftsum(1)),col(1));hold on;
        plot(x2(:,1),slip(i,1 + nftsum(1):nftsum(2)),col(2)); hold on;
        plot(x3(:,1),slip(i,1 + nftsum(2):nftsum(3)),col(3)); hold on;

        H(1)=area(x1(:,1),slip(i,1:nftsum(1)));set(H(1),'FaceColor',col(1));alpha(.3);hold on;
        H(2)=area(x2(:,1),slip(i,1 + nftsum(1):nftsum(2)));set(H(2),'FaceColor',col(2));alpha(.3);hold on;
        H(3)=area(x3(:,1),slip(i,1 + nftsum(2):nftsum(3)));set(H(3),'FaceColor',col(3));alpha(.3);hold on;
        xlim([-300 220]);

        yyaxis right;
        plot(x1(:,1),rupt(i,1:nftsum(1)),col(1), 'Linewidth', 2); hold on;
        plot(x2(:,1),rupt(i,1 + nftsum(1):nftsum(2)),col(2),'Linewidth', 2); hold on;  
        plot(x3(:,1),rupt(i,1 + nftsum(2):nftsum(3)),col(3),'Linewidth', 2); hold on;  
        ylim ([0 100]);
        
        if i == 2704 
            yyaxis left;
            ylabel('Slip (m)');
            yyaxis right;
            ylabel('Time (s)');
        end
        set(gca, 'fontweight','bold');
        set(gca,'xtick',[]);
        
        %Calculate magnitude
        
        mom = 0; len = 0;
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
        mom = alp * mom * min(len*1e3,lockdepth)*1e3;
        mag(mtag) = 2/3*log10(mom)-6.07; 
    end    
    
    if i == 2704
        subplot(nfig,1,9)
        plot(x1(1:nft(1),1),-ns(i,1:nftsum(1))*mius/1e6,col(1), 'LineWidth', 2); hold on;
        plot(x2(1:nft(2),1),-ns(i,1 + nftsum(1):nftsum(2))*mius/1e6,col(2), 'LineWidth', 2); hold on;
        plot(x3(1:nft(3),1),-ns(i,1 + nftsum(2):nftsum(3))*mius/1e6,col(3), 'LineWidth', 2); hold on;
        xlim([-300 220]);
        ylabel('Strength (MPa)');
        set(gca, 'fontweight','bold');
        set(gca,'xtick',[]);
    end
end

subplot(10,1,10)
plot(x1(:,1), x1(:,2), col(1), 'LineWidth', 2); hold on;
plot(x2(:,1), x2(:,2), col(2), 'LineWidth', 2); hold on;
plot(x3(:,1), x3(:,2), col(3), 'LineWidth', 2); hold on;
ylim([-20 70]); xlim([-300 220]);xlabel('NW-SE (km)');ylabel('SW-NE (km)');
set(gca, 'fontweight','bold');
set(gcf, 'color', 'white');

char = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)' ,'(h)', '(i)', '(j)'];
for i = 1:10
    h1(i) = subplot(nfig,1,i);
    text(h1(i),0.025,0.85,char((i-1)*3+1:i*3),'Units','normalized','FontSize',12)
end

for i = 1:8
    h1(i) = subplot(nfig,1,i);
    text(h1(i),0.025,0.55,num2str(mag(i),'%3.2f'),'Units','normalized','FontSize',12)
end

htmp = subplot(nfig,1,4);
text(htmp,0.85,0.85,'1857-like','Units','normalized','FontSize',12);
htmp = subplot(nfig,1,5);
text(htmp,0.85,0.85,'1857-like','Units','normalized','FontSize',12);
htmp = subplot(nfig,1,6);
text(htmp,0.85,0.85,'1812-like','Units','normalized','FontSize',12);
htmp = subplot(nfig,1,8);
text(htmp,0.85,0.85,'1918-like','Units','normalized','FontSize',12);
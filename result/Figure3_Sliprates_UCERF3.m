clear all; close all;
% Last modified on 10/13/2021. 
% The script is to plot the lastest version of Figure 3 with geologic slip rates.
% Dependency:
% -it calls the function observed_sliprates.m
% -./mesh/ should be available.

vert2 = load(strcat('mesh/vert.txt')); vert2 = vert2/1e3;
nsmpnv2 = load(strcat('mesh/nsmpnv.txt'));
nsmp2 = load(strcat('mesh/nsmp.txt'));

[paleoloc,rates,loc_eqdyna] = observed_sliprates;


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


for i = 1: size(x3,1)
    if x3(i,1)<(loc_eqdyna(2,1) - loc_eqdyna(1,1))/2+loc_eqdyna(1,1)
        geo_rate(i,1:3) = rates(1,1:3);
    elseif x3(i,1)>(loc_eqdyna(2,1) - loc_eqdyna(1,1))/2+loc_eqdyna(1,1) && x3(i,1)<(loc_eqdyna(3,1) - loc_eqdyna(2,1))/2+loc_eqdyna(2,1)
        geo_rate(i,1:3) = rates(2,1:3);
    else
        geo_rate(i,1:3) = rates(3,1:3);
    end
end
for i = 1:size(x1,1)
    geo1_rate(i,1:3) = rates(4,1:3);
end
for i = 1:size(x2,1)
    geo2_rate(i,1:3) = rates(5,1:3);
end

fig1=figure(1); 
set(fig1,'position', [0.1 0.1 7 5]*96);
subplot('position', [0.1 0.45 0.8 0.45]);
plot(x3(:,1), geo_rate(:,1), 'k:'); hold on;
plot(x3(:,1), geo_rate(:,2), 'r:'); hold on;
plot(x3(:,1), geo_rate(:,3), 'b:'); hold on;

plot(x1(:,1), geo1_rate(:,1), 'k:'); hold on;
plot(x1(:,1), geo1_rate(:,2), 'r:'); hold on;
plot(x1(:,1), geo1_rate(:,3), 'b:'); hold on;

plot(x2(:,1), geo2_rate(:,1), 'k:'); hold on;
plot(x2(:,1), geo2_rate(:,2), 'r:'); hold on;
plot(x2(:,1), geo2_rate(:,3), 'b:'); hold on;

for k = 1:3
    model = k;
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
    elseif model == 4
        path = './work_vis7_fs0.5_D1/'
        k = 1;
        filetag = [1,1610;];
        filenum = 2;
        mius = 0.5;
        cycstart = 1;  
        Mname = '7fs05D1';        
    elseif model == 5
        k = 1;
        filetag = [1,1256;];
        filenum = 2;
        path = './work_vis7_fs0.5/'
        mius = 0.5;
        cycstart = 1;  
        Mname = '7fs05'     
    elseif model == 6
         k = 1;
        filetag = [1,1765;];
        filenum = 2;       
        path = './work_vis7_fs0.5_lowsd/'
        mius = 0.7;
        cycstart = 1;  
        Mname = '7fs05lsd'        
    elseif model == 7
         k = 1;
        filetag = [1,1749;];
        filenum = 2;       
        path = './work_vis4.2_fs0.3_lowsd/'
        mius = 0.3;
        cycstart = 1;  
        Mname = '42fs03lsd'         
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
        ictag = ictag + 1
        ntag = ntag + 1
        i
        tmp = res((ictag-1)*ntotnd+1:ictag*ntotnd, :);
        ns(i,:) = tmp(1:ntotnd, 2);
        ss(i,:) = tmp(1:ntotnd, 1);
        slip(i,:) = tmp(1:ntotnd, 3);
        if ntag >1
            accumulateslip(ntag,:) = accumulateslip(ntag-1,:) + slip(i,:);
        end
        rupt(i,:) = tmp(1:ntotnd, 5);
        maxslip = max(slip(i,:));
    end 
    

    
    subplot('position', [0.1 0.45 0.8 0.45])
    nt = size(accumulateslip,1);
    rate = accumulateslip(nt,:)/tspan; %mm/yr
    h(k) = plot(x1(:,1),rate(1,1:nftsum(1)),col(k), 'LineWidth', 2);hold on;
    
    plot(x2(:,1),rate(1,1 + nftsum(1):nftsum(2)),col(k),'LineWidth', 2); hold on;

    plot(x3(:,1),rate(1,1 + nftsum(2):nftsum(3)),col(k),'LineWidth', 2); hold on; 
    xlim([-300 220]);ylabel('Slip rate (mm/yr)');
    set(gca,'xtick',[]);
    set(gca, 'fontsize', 12, 'fontweight', 'bold');
    %text(-40, rate(1,300)-0.5, Mname, 'FontSize',12, 'FontWeight', 'bold'); 

    subplot('position', [0.1 0.2 0.8 0.2])
    plot(x1(:,1), x1(:,2), col(1), 'LineWidth', 2); hold on;
    plot(x2(:,1), x2(:,2), col(2), 'LineWidth', 2); hold on;
    plot(x3(:,1), x3(:,2), col(3), 'LineWidth', 2); hold on;
    ylim([-20 70]); xlim([-300 220]);xlabel('NW-SE (km)');ylabel('SW-NE (km)');axis equal;
 	set(gca, 'fontsize', 12, 'fontweight', 'bold');
end

set(gcf, 'color', 'white');

%
% legend([h(1) h(2) h(3) h(4) h(5) h(6)], {'84fs03', '84fs05', '84fs07', 'uni84fs03', '126fs05', '168fs07'});
legend([h(1) h(2) h(3)], {'Model B', 'Model A', 'Model C'});
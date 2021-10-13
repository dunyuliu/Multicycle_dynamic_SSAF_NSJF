clear all; close all;
% Last modified on 10/13/2021.
% The script is to plot panels of Figure 4, 7 and 8.
% The panels show slip distributions in Model A, B and C for several kyrs.
% Parameters:
% -model = 1, 2, 3 which stand for Model B, A, C. 
% -thres defines the slip threshold. If an earthquake has maximum slip less
% than thres, it will not be shown.
% -scale is the scaling factor to show the results. Larger value will show
% higher area for slip distribution.
% -tstart defines the kyrs that you want to start the plotting.
%    In the manuscript, for Model A, periods start from 3, 6, 9, 12 kyrs.
%    For Model B, periods start from 4 and 6 kyrs.
%    For Model C, periods start from 6 and 9 kyrs. 
% Dependency:
% -mesh.mat, which is the saved mesh Matlab function space or ./mesh/.

% Current parameters are for Figure 4c. 
model = 2;
thres = 1;
scale = 30;
tstart = 9; 
fig1=figure(1); 
set(fig1,'position', [0.3 0.3 7 7]*96);

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

    k = 1;
    clear accumulateslip ttot ic tinte res;
    %%
    if model == 1 
        filenum = 3;
        filetag = [1,154,1378;];
        tend = tstart + 2; % 2kyrs interval.
        ymin = tstart - 0.5; ymax = tend + 0.2;
        path = './work_vis4.2_fs0.3/'
        mius = 0.3;
        cycstart = 1;%The first rupture of most of the system.
        Mname = '42fs03' %Model B
    elseif model == 2
        filenum = 3;
        filetag = [1,1256,2929];
        tend = tstart + 3;
        ymin = tstart - 0.5; ymax = tend + 0.2;       
        path = './work_vis7_fs0.5/'
        mius = 0.5;
        cycstart = 1;
        Mname = '7fs05' %Model A
    elseif model == 3
        filenum = 3;
        filetag = [1,1436,3284];
        tend = tstart + 3;
        ymin = tstart - 0.5; ymax = tend + 0.2;         
        path = './work_vis12_fs0.7/'
        mius = 0.7;
        cycstart = 1;  
        Mname = '12fs07' %Model C        
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
    col = ['b','r','k','m','g','c','y'];
    col2 = ['b:','r:','k:','m:','g','c','y'];

    for kk = 1: filenum
        if kk == 1
            ic = load(strcat(path,'cyclelog.txt',num2str(filetag(k,kk))));
            res = load(strcat(path,'totalop.txt',num2str(filetag(k,kk))));
            tinte  = load(strcat(path,'interval.txt',num2str(filetag(k,kk))));
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
    
    ntag = 0;
    ictag = 0;
    tspan = ttot(ic(2)-ic(1)+1) - ttot(1);    
    for i = ic(1): ic(2)
        ictag = ictag + 1;
        if ttot(i)>=tstart && ttot(i)<=tend
            ntag = ntag + 1
            if (mod(i,100) == 1)
                i
            end
            tmp = res((ictag-1)*ntotnd+1:ictag*ntotnd, :);
            ns(i,:) = tmp(1:ntotnd, 2);
            ss(i,:) = tmp(1:ntotnd, 1);
            slip(i,:) = tmp(1:ntotnd, 3);
            maxslip = max(slip(i,:));

            if (maxslip > thres) 
                xx1 = [x1(:,1)',fliplr(x1(:,1)')]; yy1 = [slip(i,1:nftsum(1))/scale + ttot(i),fliplr(zeros(1,nftsum(1)) + ttot(i))];
                H(1)= fill(xx1,yy1,col(1)); alpha(.3);hold on;
                xx2 = [x2(:,1)',fliplr(x2(:,1)')]; yy2 = [slip(i,1 + nftsum(1):nftsum(2))/scale + ttot(i),fliplr(zeros(1,nft(2)) + ttot(i))];
                H(2)= fill(xx2,yy2,col(2)); alpha(.3);hold on; 
                xx3 = [x3(:,1)',fliplr(x3(:,1)')]; yy3 = [slip(i,1 + nftsum(2):nftsum(3))/scale + ttot(i),fliplr(zeros(1,nft(3)) + ttot(i))];
                H(3)= fill(xx3,yy3,col(3)); alpha(.3);hold on;  
                text(152 + mod(ntag,3)*15, ttot(i), strcat(num2str(i)), 'FontSize', 6, 'Fontweight','bold');
            end
        end
    end
    plot(x1(:,1), x1(:,2)/150 + tstart-0.4, col(1), 'LineWidth', 2); hold on;
    plot(x2(:,1), x2(:,2)/150 + tstart-0.4, col(2), 'LineWidth', 2); hold on;
    plot(x3(:,1), x3(:,2)/150 + tstart-0.4, col(3), 'LineWidth', 2); hold on;
    
    xlim([-250 160]);ylim([ymin ymax]);ylabel('Time (kyrs)');
    xlabel('NW-SE (km)');
    
    xxx=[120, 120]; yyy = [ymin+0.05,ymin+0.05+5/scale];
    line(xxx,yyy,'linewidth',3,'color','k');
    text(125,ymin+0.05+5/scale/2,'slip = 5 m');
    
set(gcf, 'color', 'white');
set(gca, 'fontsize', 12, 'fontweight', 'bold');

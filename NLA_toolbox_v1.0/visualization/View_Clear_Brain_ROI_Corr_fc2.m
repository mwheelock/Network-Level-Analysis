function pairsOut=View_Clear_Brain_ROI_Corr_fc2(Anat,IM,M,dataOut,fc,Pmax,params)
%
% This function draws ROIs on the brain and connecting rods colored by sign
% of correlation values in C. This can be the fc correlation itself if that 
% is passed, or it can be the correlation of brain fc with behavior if that 
% matrix is passed. Pairs of ROIs to be included are taken from the
% elements below threshold (Pmax) in matrix Proi.
% Inputs:   Anat - anatomy for visualization of hemispheres
%           IM   - ROI-Network info, including ROIxyz field
%           M    - either the matrix of network pair tests, or the
%                   threshold for the network-pair p-values, or a list of
%                   Nx2 Network pairs
%           Pnn  - p-val matrix for thresholding network pairs
%           Proi - p-val matrix for thresholding ROI pairs
%           fc    - correlation matrix (fc or fcBx)
%           Pmax - is the Proi threshold for significance.
%           toDraw - flag for drawing ROIs (1) or just outputting roi pairs

%% Parameters
params.brain=1;
params.fig=0;
if ~isfield(params,'toDraw'), params.toDraw=1;end
if ~isfield(params,'ScaleRadius'), params.ScaleRadius=0;end
if ~isfield(params,'type'), params.type='perc';end
if (strcmp(params.type,'Spearman')||strcmp(params.type,'Pearson'))
    params.type='perc';
end
FC.tp1=fc{1};
FC.tp2=fc{2};
clear fc
cMapN=hsv(2800);cMapN=cMapN(1401:2400,:);
cMapP=cat(1,summer(500),flip(autumn(500),1));
    
%% Select Networks of interest
NNsig=M;
NnnSig=size(NNsig,1);
pairsOut={};


%% Select Pairs of interest
if NnnSig>0

for j=1:NnnSig
    
    mask=zeros(size(IM.key,1));
    mask(IM.key(:,2)==NNsig(j,1),IM.key(:,2)==NNsig(j,2))=1;
    % Grab sig ROI pairs for each time step
    for ts=1:2
    pairsOut{j,ts}=[];    
    Proi=squeeze(dataOut.pval(:,:,ts));
    [pairsOut{j,ts}(:,1),pairsOut{j,ts}(:,2)]=...
        ind2sub(size(Proi),find(tril(mask,-1).*(Proi<=Pmax)));
%     pairsOut{j,ts}=pairs;
    if ~isempty(pairsOut)
    
Npairs=size(pairsOut{j,ts},1);
Rs=unique(pairsOut{j,ts}(:));
Nrois=size(Rs,1);

    
radius=params.roiradius;
roi.radius=repmat(radius,Nrois,1); %%## Set radius (volume?) to proportional
roi.coord=zeros(Nrois,3);
roi.color=zeros(Nrois,3);
Conn=zeros(Npairs,6);
    
% Set parameters for ROI locations and colors
for k=1:Nrois
    roi.coord(k,:)=IM.ROIxyz(Rs(k),:);
    roi.color(k,:)=IM.cMap(IM.key(Rs(k),2),:);
    if params.ScaleRadius
        temp = radius+sum(pairsOut{j,ts}(:)==Rs(k));
        if temp<=10;
    roi.radius(k,1)=radius+sum(pairsOut{j,ts}(:)==Rs(k)); % Scale Rad, Area, Vol?
        else
            roi.radius(k,1)=10;
        end
    end
end

% Set parameters for connecting sticks
for k=1:Npairs
    Conn(k,:)=[find(Rs==pairsOut{j,ts}(k,1)),...
        find(Rs==pairsOut{j,ts}(k,2)),0,0,0,0];
    rho=dataOut.rho(Rs(Conn(k,1)),Rs(Conn(k,2)),ts)>0;
    fc=squeeze(FC.(['tp',num2str(ts)])(Rs(Conn(k,1)),Rs(Conn(k,2)),:));
    switch params.type
        case 'perc'
    p=round(1000*(sum(fc>0)/length(fc)));if p==0,p=1;end
    pairsOut{j,ts}(k,4)=sum(fc>0)/length(fc);
        case 'mean'
    p=round(1000*(mean(fc))+500);if p<=0,p=1;elseif p>1000,p=1000;end
    pairsOut{j,ts}(k,4)=mean(fc);
        case 'median'
    p=round(1000*(median(fc))+500);if p<=0,p=1;elseif p>1000,p=1000;end
    pairsOut{j,ts}(k,4)=median(fc);
    end
    pairsOut{j,ts}(k,3)=dataOut.rho(Rs(Conn(k,1)),Rs(Conn(k,2)),ts);
     pairsOut{j,ts}(k,5)=dataOut.pval(Rs(Conn(k,1)),Rs(Conn(k,2)),ts);
    pairsOut{j,ts}(k,6)=min(fc);
    pairsOut{j,ts}(k,7)=max(fc);
    pairsOut{j,ts}(k,8)=(mean(fc)/std(fc))*sqrt(length(fc)-1);
    pairsOut{j,ts}(k,9)=mean(fc);
    pairsOut{j,ts}(k,10)=std(fc);
    pairsOut{j,ts}(k,11)=median(fc);
    Ns=length(fc);
    sfc=sort(fc);
    pairsOut{j,ts}(k,12)=sfc(round(Ns/4));
    pairsOut{j,ts}(k,13)=sfc(round(3*(Ns/4)));
    
    % Set colors of rod
    switch rho
        case 1  % positive fcBx correlation
            Conn(k,3:5)=cMapP(p,:);
        case 0  % negative fcBx correlation
            Conn(k,3:5)=cMapN(p,:);
    end
end


% Draw cool stuff on brains
if params.toDraw   
    if size(pairsOut{ts},1)>=1
figure('Color','w');fcBx_PercFC_Scat2(pairsOut{ts},IM); %perc f>0
end
disp(['Drawing ROIs for ',IM.Nets{NNsig(j,1)},...
    ' and ',IM.Nets{NNsig(j,2)}])
Draw_ROIs_Through_Cortex_3_Views(Anat,roi,Conn,params);
%end
ax = findall(gcf, 'Type', 'axes');
subplot(ax(2))
title([{params.group{ts}};...
    [{[IM.Nets{NNsig(j,1)},'-',IM.Nets{NNsig(j,2)},...
    ', p < ',num2str(Pmax)]}]])

    else
        disp(['No ROI pairs with p < ',num2str(Pmax)])
    end
    end
    end
end
end    
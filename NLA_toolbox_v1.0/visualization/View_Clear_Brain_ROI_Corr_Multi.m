function View_Clear_Brain_ROI_Corr_Multi(Anat,IM,M,Pnn,Proi,C,Pmax,params,g)
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
%           C    - correlation matrix (fc or fcBx)
%           Pmax - is the Proi threshold for significance.
%           toDraw - flag for drawing ROIs (1) or just outputting roi pairs

%% Parameters
params.brain=1;
params.fig=0;
radius = params.roiradius;
if ~isfield(params,'toDraw'), params.toDraw=1;end
if ~isfield(params,'ScaleRadius'), params.ScaleRadius=0;end

    
%% Select Networks of interest
NNsig=[];
if size(M,2)>2         % M is the matrix of desired network pairs
[NNsig(:,1),NNsig(:,2)]=ind2sub(size(M),find(M==1));

elseif numel(M)==1   % M is the threshold
[NNsig(:,1),NNsig(:,2)]=ind2sub(size(Pnn),find(tril(Pnn<=M,-1)));

else                   % M is an Nx2 array of Network pairs
     NNsig=M;
end
NnnSig=size(NNsig,1);


%% Select Pairs of interest
disp(['Drawing up to ',num2str(NnnSig),' sets of brains'])

pairs=[];
n=0;
for j=1:NnnSig
    temp=[];
    mask=zeros(size(IM.key,1));
    mask(IM.key(:,2)==NNsig(j,1),IM.key(:,2)==NNsig(j,2))=1;
    [temp(:,1),temp(:,2)]=...
        ind2sub(size(Proi),find(tril(mask,-1).*(Proi<=Pmax)));
    pairs=cat(1,pairs,temp);
end    
Npairs=size(pairs,1);
Rs=unique(pairs(:));
Nrois=size(Rs,1);
roi.radius=repmat(radius,Nrois,1);
if ~isfield(roi,'radius'), roi.radius=repmat(5,Nrois,1); %%## Set radius (volume?) to proportional
end
roi.coord=zeros(Nrois,3);
roi.color=zeros(Nrois,3);
Conn=zeros(Npairs,3);
    
% Set parameters for ROI locations and colors
for k=1:Nrois
    roi.coord(k,:)=IM.ROIxyz(Rs(k),:);
    roi.color(k,:)=IM.cMap(IM.key(Rs(k),2),:);
    if params.ScaleRadius
    roi.radius(k,1)=4+sum(pairs(:)==Rs(k)); % Scale Rad, Area, Vol?
    end
end

% Set parameters for connecting sticks
for k=1:Npairs
    Conn(k,:)=[find(Rs==pairs(k,1)),find(Rs==pairs(k,2)),...
        C(pairs(k,1),pairs(k,2))];
end

% Draw cool stuff on brains
Draw_ROIs_Through_Cortex_3_Views(Anat,roi,Conn,params); 
ax = findall(gcf, 'Type', 'axes');
subplot(ax(2))
title([{params.group{g}};...
    [{[IM.Nets{NNsig(j,1)},'-',IM.Nets{NNsig(j,2)},...
    ', p < ',num2str(Pmax)]}]])
   
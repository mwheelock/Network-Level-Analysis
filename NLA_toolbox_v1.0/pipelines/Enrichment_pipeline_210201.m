% Script for running Enrichment analyses between fc data and behavioral
% data.
% It is assumed that the fc data is organized as ROI-ROI-Subject and the
% behavioral data is a vector with the same order and number of subjects
% in the fc array.

outputdir = 'filepath to output directory/';
addpath(genpath('filepath to toolbox and subfolders'));

%% Load in data

% Load behavioral data
Bxname='behavior name';
group = 'group name';
Behavior=importdata(strcat(outputdir,'filepath to behavioral data'));

fc = load('name of fc data');
Bx = 'behavior name';

%Load fc data
corrmat=load('filepath to your fc data')

% Load IM structure containing model for fc structure
load('filepath to IM structure');

%% plot average fc organized by networks

% if z-values
    rmatAve=FisherZ2R(mean(FisherR2Z(corrmat),3)); 
% if r-values
     rmatAve=mean(corrmat,3);

% % load anat (surface file) % %
load('MNI152nl_on_TT_coord_meshes_32k','MNIl','MNIr'); % MNI152 template in TT
%load('MNI_coord_meshes_32k','MNIl','MNIr'); % MNI152 in MNI
Anat.CtxL=MNIl;Anat.CtxR=MNIr;
Anat.alpha = 0.25; % how see through is mesh?
clear MNIl MNIr

% % Visualize rmat organized by IM file % %
Matrix_Org3(rmatAve(IM.order,IM.order),...
    IM.key,10,[-0.3,0.3],IM.cMap,1);

% % Visualize IM Model on Brain with Network Names and Colors 
params.radius = 4; %specify ROI size for display on brain surface
Anat.ctx='std'; %typically use 'inf' for adult mesh, mesh with CB is 'std'
View_ROI_Modules(IM,Anat,IM.ROIxyz,params);

% % plot ROI count in each network % %
figure;histogram(IM.key(:,2));title(strrep(IM.name,'_',' '))
set(gca,'XTick',[1:max(IM.key(:,2))],'XTickLabel',IM.Nets);ylabel('Nrois')

%% Fishers z transform ?

%if r values
fc=FisherR2Z(corrmat);
%if z values
fc=corrmat;

%% Set parameters
params.np=1e4;        % Number of permutations: set to 1 to get simple stats
                      % Set to 1e4 to get reasonable distribution.
params.type='Pearson'; % flavor of correlations: 'Pearson' or 'Spearman'
params.Pmax=0.01;
params.nnPmax=0.05; 
params.group = [group {}];
params.BxName=Bxname;
params.B=1; % Bonferonni level (i.e., how many Bx are you using?)
params.fn=strcat(Bxname,'_',num2str(params.np),'_',num2str(params.Pmax),'_',IM.name,'.mat'); 

%% Run Enrichment analyses

% assumes fc are Fisher z transformed r values
dataOut=fcBx_Enrich_1tp(fc,Bx,IM,params);
cd(outputdir)
save(params.fn,'dataOut','fc','Bx','params');
clear dataOut

params.B = 2;
% To Re-visualize data after loading dataOut structure 
Perm_1tp_Figs2(dataOut,IM,params.B,params);

%% make stats tables

RC = [5,4;10,2;8,2;12,5]; % specify row,column of network statistic to extract

Nnets=length(IM.Nets);
sizenets=dataOut.Np;
for j=1:length(dataOut.TNidx);
    nsize(dataOut.TNidx(j))=sizenets(j);
end
netsize=reshape(nsize,Nnets,Nnets,[]);

params.roiradius=5;
params.toDraw=0;
for j=1:length(RC)    
pairsOut{j}=View_Clear_Brain_ROI_Corr_fc2(Anat,dataIN.IM,RC(j,:),...
    dataOut,dataIN.fc,0.05,params);
end
for j=1:size(RC,1)
    for g=1:length(params.group);
        name=params.group(g);name=name{1,1};
NetStats.nets{j} = strcat(dataIN.IM.Nets{RC(j,1)},'-',dataIN.IM.Nets{RC(j,2)});
NetStats.(name).Chipvals(j) =dataOut.Chi_EWpval(RC(j,1),RC(j,2),g);
NetStats.(name).Chistats(j) =dataOut.Chi_stats(RC(j,1),RC(j,2),g); 
NetStats.(name).Nethits(j)=size(pairsOut{1,j}{1,1},g);
NetStats.(name).N(j) = netsize(RC(j,1),RC(j,2));
NetStats.(name).phieffectsize(j)=sqrt(dataOut.Chi_stats(RC(j,1),RC(j,2),g)/(netsize(RC(j,1),RC(j,2))));
    end
end
save NetworkStatistics.mat NetStats

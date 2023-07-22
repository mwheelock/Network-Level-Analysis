% Example 4: get parcel files for visualization
% assumes that the cifti files is cortex only with 59412 vertices (fsLR32k)
addpath(genpath(pwd));
addpath(genpath('/data/wheelock/data1/software/cifti-matlab-master')); % downloaded from https://www.mathworks.com/matlabcentral/fileexchange/56783-washington-university-cifti-matlab

outdir = '/data/wheelock/data1/people/Cindy/BCP/ParcelPlots' % Change me % 

parcelname = 'Parcellation1'% Change me %
parcels_path = '/data/wheelock/data1/parcellations/DustinScheinostNeonateParcellation/baby_atlas_in_MNI_space_fs_LR32k_dilated10mm_cleanedpt1thresh_LR_minsize10.dlabel.nii'; % Change me % 

parceldata = cifti_read(parcels_path);
Linds=with_without_mw_conversion('Lindfull');
Rinds=with_without_mw_conversion('Rindfull');

n_verts_per_hem = 32492;
Parcels.CtxL = zeros(n_verts_per_hem,1);% total vertices
Parcels.CtxR = zeros(n_verts_per_hem,1);

Parcels.CtxL(Linds) = parceldata.cdata(1:length(Linds));
Parcels.CtxR(Rinds-n_verts_per_hem) = parceldata.cdata(length(Linds)+1:end);

% matpath= fullfile(outdir,strcat('/Parcels_',parcelname,'.mat'));
% save(matpath,'Parcels');

%% View all parcels on MNI
load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr

Nparcels = length(setdiff(unique([Parcels.CtxL;Parcels.CtxR]),0))

% (1) Add parcel nodes to Cortices
% Parcels.CtxL(Parcels.CtxL==0) = Nparcels+1; % set the last one to black to draw borders (if available) and medial wall
% Parcels.CtxR(Parcels.CtxR==0) = Nparcels+1;
Anat.CtxL.data=Parcels.CtxL;
Anat.CtxR.data=Parcels.CtxR;
% (2) Set parameters to view as desired
params.Cmap.P=colorcube(Nparcels);
params.Cmap.P(Nparcels+1,:) = [0 0 0];
params.TC=1;
params.ctx='inf';           % 'std','inf','vinf'
figure;
ax = subplot(2,1,1);
params.fig_handle = ax;
params.view= 'lat';       % 'dorsal','post','lat','med'
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);
title(parcelname,'interpreter','none','color','k')
ax = subplot(2,1,2);
params.fig_handle = ax;
params.view ='med';
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);

set(gcf,'color','w','InvertHardCopy','off');
% print('Example3.png','-dpng');
% print(gcf,fullfile(outdir,[parcelname,'.png']),'-dpng');
return


%% (optional) obtain distance matrices
neighbors_path = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_surf_neighbors_LR_normalwall.mat'; %Change me%
load(neighbors_path,'neighbors');

dist_path = {'/data/wheelock/data1/parcellations/SurfaceFiles/distmat_surf_geodesic_L_uint8.mat',...
    '/data/wheelock/data1/parcellations/SurfaceFiles/distmat_surf_geodesic_R_uint8.mat'};%Change me%

dist(1) = load(dist_path{1});
dist(2) = load(dist_path{2});
distL = squareform(dist(1).D);
distR = squareform(dist(2).D);

fcn_calc_neighs_per_parcel(parcels_path,distL,distR,neighbors,matpath)


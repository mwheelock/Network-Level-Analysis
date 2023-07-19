% N.B. one problem is that the Gordon parcellation and Myers parcellation
% are not covering the surface
% load('/data/nil-bluearc/GMT/Evan/Atlases/32k_ConteAtlas_v2_distribute/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat');
% distL = distances(Lindtrunc,Lindtrunc);
% distR = distances(Rindtrunc,Rindtrunc);
function fcn_calc_neighs_per_parcel(parcels_path,distL,distR,neighbors,matpath)
Lindtrunc = with_without_mw_conversion('Lindtrunc');
Rindtrunc = with_without_mw_conversion('Rindtrunc');
Lindfull = with_without_mw_conversion('Lindfull');
Rindfull = with_without_mw_conversion('Rindfull');
assert(length(Lindtrunc)==length(distL));
assert(length(Rindtrunc)==length(distR))
%% Load parcellation
parcel = cifti_read(parcels_path);
parcel = parcel.cdata;
parcel_ids = setdiff(unique(parcel),0);
parcelL = parcel(Lindtrunc);
parcelR = parcel(max(Lindtrunc)+1:end);
parcelL_ids = setdiff(unique(parcelL),0);
parcelR_ids = setdiff(unique(parcelR),0);
%% Load vertices xyz
% load('/data/wheelock/data1/people/Cindy/BCP/ParcelPlots/MNI_coord_meshes_32k.mat');
% Lnodes = MNIl.nodes;
% Rnodes = MNIr.nodes;
% Lnodes = gifti('/data/wheelock/data1/parcellations/Arslan2018/Scripts/surface/Conte69.L.midthickness.32k_fs_LR.surf.gii');
Lnodes = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii');
Lnodes = Lnodes.vertices;
% Rnodes = gifti('/data/wheelock/data1/parcellations/Arslan2018/Scripts/surface/Conte69.R.midthickness.32k_fs_LR.surf.gii');
Rnodes = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.midthickness.32k_fs_LR.surf.gii');
Rnodes = Rnodes.vertices;

%% distance between centroid of parcels
disp('Calculating centroid distances')
all_parcel_centroids = NaN(size(parcel_ids));
for i = parcelL_ids'
    parcel_verts = find(parcelL==i);
    this_parcel_dmat = distL(parcel_verts,parcel_verts);
    summed_distances = sum(this_parcel_dmat,2);
    [minsummeddist, centroidindex] = min(summed_distances);
    all_parcel_centroids(i) =  parcel_verts(centroidindex);
end
for i = parcelR_ids'
    parcel_verts = find(parcelR==i);
    this_parcel_dmat = distR(parcel_verts,parcel_verts);
    summed_distances = sum(this_parcel_dmat,2);
    [minsummeddist, centroidindex] = min(summed_distances);
    all_parcel_centroids(i) =  parcel_verts(centroidindex);
end

parcels_dmat = uint8(ones(length(parcel_ids))*255);
parcels_dmat(parcelL_ids,parcelL_ids) = distL(all_parcel_centroids(parcelL_ids),all_parcel_centroids(parcelL_ids));
parcels_dmat(parcelR_ids,parcelR_ids) = distR(all_parcel_centroids(parcelR_ids),all_parcel_centroids(parcelR_ids));

ROIxyz = single(NaN(length(parcel_ids),3));
ROIxyz(parcelL_ids,:) = Lnodes(Lindfull(all_parcel_centroids(parcelL_ids)),:);
ROIxyz(parcelR_ids,:) = Rnodes(Rindfull(all_parcel_centroids(parcelR_ids))-length(Lnodes),:);

% save(['/data/wheelock/data1/people/Cindy/BCP/ParcelPlots/Parcels_',parcel_name,'.mat'],'parcels_dmat','ROIxyz','-append')
% return


%% calculate neighs for parcel
disp('Calculating neighbors')
neighs = uint8(NaN(max(parcel)));
for i = setdiff(unique(parcel)',0)
    idx = neighbors(parcel==i,2:end);
    idx = unique(idx(~isnan(idx)));
    parcelneighs = setdiff(unique(parcel(idx)),[0,i]);
    if ~isempty(parcelneighs)
        neighs(i,parcelneighs) = 1;
    end
end
%% calculate minimum distance for parcel
disp('Calculating minimum distances')
mindist = uint8(NaN(max(parcel)));
for i = setdiff(unique(parcelL),0)'
    for j = i+1:max(parcelL)
        tmp = distL(parcelL==i,parcelL==j);
        mindist(i,j) = min(tmp(:));
    end
end
for i = setdiff(unique(parcelR),0)'
    for j = i+1:max(parcelR)
        tmp = distR(parcelR==i,parcelR==j);
        mindist(i,j) = min(tmp(:));
    end
end
mindist = mindist+mindist';
mindist(max(unique(parcelL))+1:end,1:max(unique(parcelL))) = 255;
mindist(1:max(unique(parcelL)),max(unique(parcelL))+1:end) = 255;
%%
save(matpath,'mindist','neighs','parcels_dmat','ROIxyz','-append')
end


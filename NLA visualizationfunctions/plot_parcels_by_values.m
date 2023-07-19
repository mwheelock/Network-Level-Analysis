% This helps to visualize homogeneity/variance at individual parcel
function plot_parcels_by_values(x,Anat,parcelview,Parcels,clim,cmap)
    %% Find parcels to remove
    idx = find(isnan(x));
    Parcels.CtxL(any(Parcels.CtxL==idx',2))=0;
    Parcels.CtxR(any(Parcels.CtxR==idx',2))=0;
    %% View all parcels on mesh
    % (1) Add parcel nodes to Cortices
    Anat.CtxL.data=Parcels.CtxL;
    Anat.CtxR.data=Parcels.CtxR;
    % (2) Set parameters to view as desired
    params.Cmap.P = NaN(length(x),3); % can be NaN too just leave these indices as gaps %I tried applycmap for specific maps but does not seem to work?
    params.Cmap.P(~isnan(x),:)=value_to_cmap(x(~isnan(x)),clim(1),clim(2),cmap);
    params.TC= 1;
    params.ctx='inf';           % 'std','inf','vinf'
    params.view= parcelview;       % 'dorsal','post','lat','med'
    params.fig_handle = gca;
    PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);
    set(gcf,'color','w');
end

function xcolors = value_to_cmap(x,cmin,cmax,cmap)
    maptounity = @(v)(v-cmin)/(cmax-cmin);
    x = round(maptounity(x)*100)+1;
    x(x>100) = 100;
    x(x<1)=1;
    xcolors = cmap(x,:);
end
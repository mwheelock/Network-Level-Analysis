function Draw_ROIs_on_Cortex(Anat,ROI)
%
% This function will draw ROI spheres on a model of the cortical surface.
% Anat is a structure that describes the anatomy.
% ROI is a structure that describes the regions of interest. 
%
% Anat fields:
%   CtxL    mesh for left cortex. mesh includes 2 fields, nodes and
%           elements.
%   CtxR    mesh for right cortex
%   view    perspective of image {'lat'(default),'med','post','dorsal',)}
%          (optional)
%   ctx    cortical surface type. {'std'(default),'inf','vinf'} (optional)
% ROI fields:
%   coord   [x,y,z] coordinates for the ROIs. It is assumed that the ROIs
%           are in the same coordinate space as the anatomy.
%   radius  size of spheres to be drawn (default = 5). (optional)
%   color   color of each sphere. Typically this matches some Network
%           organization.  If this is also not present, all
%           spheres will be colored white. (optional)
%   Network Nx2 array of ROI indices (1st col) and network membership (2nd
%           column). (optional)


%% Under the hood Parameters
if ~isfield(Anat,'view'),Anat.view='lat';end
if ~isfield(Anat,'ctx'),Anat.ctx='std';end
Anat.CtxL.data=zeros(size(Anat.CtxL.nodes,1),1);
Anat.CtxR.data=zeros(size(Anat.CtxR.nodes,1),1);

Nroi=size(ROI.coord,1);
% radius= ROI.roiradius;
% ROI.radius = repmat(radius,[Nroi,1]);
if ~isfield(ROI,'radius'), ROI.radius=repmat(5,[Nroi,1]);end
if ~isfield(ROI,'Network')
    ROI.Network=[[1:Nroi]',ones(Nroi,1)];
    Nnet=1;
else
    Nnet=max(ROI.Network(:,2));
end
if ~isfield(ROI,'color')
    if Nnet>1
        Cmap=jet(Nnet);
        ROI.color=zeros(Nroi,3);
        for j=1:Nnet
           keep=find(ROI.Network(:,2)==j);
           N=length(keep);
           ROI.color(keep,:)=repmat(Cmap(j,:),[N,1]);
        end
    else
        ROI.color=1.0.*ones(Nroi,3);
    end
end

if ~exist('params','var')
params.Scale=1;
params.Th.P=1;
params.Th.N=-1;
params.lighting='none';
end
if isfield(Anat,'alpha'),params.alpha=Anat.alpha;end
% figCol='w'; % make fig white for easy printing. make black for good pptx


%% Draw cortex
Image_Mesh_Data_2hem(Anat.CtxL,Anat.CtxR,params,Anat.view,Anat.ctx,0);
colorbar off
hold on


%% Draw ROIs
foci.lighting=params.lighting;
for j=1:Nnet
    keep=find(ROI.Network(:,2)==j);
    foci.color=ROI.color(keep,:);
    foci.radius=ROI.radius(keep,:);
    foci.location=ROI.coord(keep,:);
    foci=AdjustFoci(foci,Anat.CtxL,Anat.CtxR,Anat.view,Anat.ctx);
    hold on
    Draw_Foci(foci,10)
end
% set(gcf,'Color',figCol)
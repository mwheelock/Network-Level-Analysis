function Draw_ROIs_Through_Cortex(Anat,ROI,Conn,params)
%
% This function will draw ROI spheres on a model of the cortical surface.
% The ROI spheres may be connected by tubes that are sized by correlation
% strength and colored by sign. The cortical surface with have a
% transparency.
% Anat is a structure that describes the anatomy.
% ROI is a structure that describes the regions of interest. 
% Conn is an Nx3 array for the connectivity to be visualized.
%
% Anat fields:
%   CtxL    mesh for left cortex. mesh includes 2 fields, nodes and
%           elements.
%   CtxR    mesh for right cortex
%   view    perspective of image {'dorsal'(default),'med','post')}
%          (optional)
%   ctx    cortical surface type. {'std'(default),'inf','vinf'} (optional)
% ROI fields:
%   coord   [x,y,z] coordinates for the ROIs. It is assumed that the ROIs
%           are in the same coordinate space as the anatomy.
%   color   color of each sphere. Typically this matches some Network
%           organization. 
%   radius  size of spheres to be drawn (default = 5). (optional)
% Conn:     columns 1&2 are ROI numbers and column 3 holds the correlation
%           value that determines the tube sizing and coloring.
%           if Conn has 5 columns rather than 3, col 3:5 are the RGB values
% params.fig    {0,1} to draw a new figure. 1 == draw new fig (default).
% params.brain  {0,1} to drawe the cortex. set to zero if figure already
%           exists and the cortex is already drawn.


%% Under the hood Parameters
if ~isfield(Anat,'view'),Anat.view='dorsal';end
if ~isfield(Anat,'ctx'),Anat.ctx='std';end
Anat.CtxL.data=zeros(size(Anat.CtxL.nodes,1),1);
Anat.CtxR.data=zeros(size(Anat.CtxR.nodes,1),1);

Nroi=size(ROI.coord,1);
if ~isfield(ROI,'radius'), ROI.radius=repmat(5,[Nroi,1]);end

Ntube=size(Conn,1);

if ~isfield(params,'fig'),params.fig=1;end
if ~isfield(params,'brain'),params.brain=1;end
if ~isfield(params,'lighting'),params.lighting='phong';end
params.Scale=1;
params.Th.P=1;
params.Th.N=-1;
params.alpha=0.25;
% params.alpha=0.025;
foci.lighting=params.lighting;
figCol='w'; % make fig white for easy printing. make black for good pptx


%% Draw cortex
if params.brain
Image_Mesh_Data_2hem(Anat.CtxL,Anat.CtxR,params,Anat.view,Anat.ctx,...
    params.fig);
colorbar off
end
hold on


%% Draw tubes
ends.location=cat(1,ROI.coord(Conn(:,1),:),ROI.coord(Conn(:,2),:));
ends=AdjustFoci(ends,Anat.CtxL,Anat.CtxR,Anat.view,Anat.ctx);
    
for j=1:Ntube
    if size(Conn,2)==3
    switch sign(Conn(j,3))
        case 1
            Tcolor=[1 0 0];
        case -1
            Tcolor=[0 0 1];
    end
    plot3(ends.location([j,j+Ntube],1),ends.location([j,j+Ntube],2),...
        ends.location([j,j+Ntube],3),...
        'Color',Tcolor,'LineWidth',5);%ceil(10*abs(Conn(j,3)))
    elseif Conn(j,6)==0
        plot3(ends.location([j,j+Ntube],1),ends.location([j,j+Ntube],2),...
        ends.location([j,j+Ntube],3),...
        'Color',squeeze(Conn(j,3:5)),'LineWidth',5);%
    elseif Conn(j,6)==1
        plot3(ends.location([j,j+Ntube],1),ends.location([j,j+Ntube],2),...
        ends.location([j,j+Ntube],3),'--',...
        'Color',squeeze(Conn(j,3:5)),'LineWidth',5);%
    end
    
end


%% Draw ROIs
foci.color=ROI.color;
foci.radius=ROI.radius;
foci.location=ROI.coord;
foci=AdjustFoci(foci,Anat.CtxL,Anat.CtxR,Anat.view,Anat.ctx);
hold on
Draw_Foci(foci,10)
set(gcf,'Color',figCol)
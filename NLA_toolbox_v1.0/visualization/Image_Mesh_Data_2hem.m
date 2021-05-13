function Image_Mesh_Data_2hem(meshL,meshR,params,side,ctx,fig)

% This function renders a mesh with the shading/coloring determined by the
% values of the data at those points.
% Update 170222, info is now needed to account for acquision direction.
% Prev, his was optimized for transverse only. now sagittal is supported.

%% Parameters
if ~isfield(params,'Scale') 
    params.Scale=0.9*max(cat(1,meshL.data(:),meshR.data(:)));end
if ~isfield(params,'Th')
    params.Th.P=0.25*params.Scale;
    params.Th.N=-params.Th.P;
end
if ~isfield(params,'Cmap'), params.Cmap='jet';end
if ~isfield(params,'OL'), params.OL=0;end
if ~isfield(params,'PD'), params.PD=0;end
if ~isfield(params,'alpha'), params.alpha=1;end % Transparency
if ~isfield(params,'TC'),params.TC=0;end  % fixed 'True-Color' colormap
if ~exist('ctx','var'), ctx='std';end
if ~exist('fig','var'), fig=1;end


%% Draw figure
if fig==1,figure('Units','inches');end
set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);


%% Image hemispheres
% Re-position meshes for standard transverse orientation etc.
[Rnodes,Lnodes]=Adjust_Brain_Pos_for_Vis(meshL,meshR,params,side,ctx);

% Set lighting and persepctive
switch side
    case {'lat','med'}
            view([-90,0]);
        light('Position',[-100,200,0],'Style','local');
        light('Position',[-50,-500,100],'Style','infinite'); % These two lines create minimal lighting good luck.
        light('Position',[-50,0,0],'Style','infinite');

if fig==1
% chbar=colorbar('location','South','Position',[0.44 0.325 0.15 0.02],...
%     'XTickLabel',{[num2str(-params.Scale)],'0.0',[num2str(params.Scale)]});
end 
    
    case {'post','dorsal'} % [x:L->R,y:P->A,z:V->D]
        if strcmp(side,'post'),view([0 0]);end
        if strcmp(side,'dorsal')
            view([0 90]);
        light('Position',[100,300,100],'Style','infinite');
        end 
        light('Position',[-500,-20,0],'Style','local');
        light('Position',[500,-20,0],'Style','local');
        light('Position',[0,-200,50],'Style','local');
end

% Image Left
dataL=CMapMix2(meshL.data,params);
FaceColor='interp';
hL=patch('Faces',meshL.elements(:,1:3),'Vertices',Lnodes,...
    'EdgeColor','none','FaceColor',FaceColor,'FaceVertexCData',dataL,...
    'FaceLighting','gouraud','FaceAlpha',params.alpha,...
    'AmbientStrength',0.25,'DiffuseStrength',.75,'SpecularStrength',.1);

% Image Right
dataR=CMapMix2(meshR.data,params);
FaceColor='interp';
hR=patch('Faces',meshR.elements(:,1:3),'Vertices',Rnodes,...
    'EdgeColor','none','FaceColor',FaceColor,'FaceVertexCData',dataR,...
    'FaceLighting','gouraud','FaceAlpha',params.alpha,...
    'AmbientStrength',0.25,'DiffuseStrength',.75,'SpecularStrength',.1);
    
axis image;
axis off
% axis vis3d
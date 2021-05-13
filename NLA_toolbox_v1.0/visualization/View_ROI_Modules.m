function View_ROI_Modules(IM,Anat,ROIxyz,params)
%
% This function generates a plot of ROI Sortings for
% various kden, scanning through the kden chosen.

%% Set up paramters
Modules=IM.key(:,2);
Cmap=IM.cMap;
ROI.Nfaces=50;
Nroi=size(Modules,1);

Params.Cmap=Cmap;
Params.Scale=length(Cmap);Params.OL=1;Params.PD=1;Params.TC=1;
Params.Th.P=0.001;Params.Th.N=-Params.Th.P;
Params.lighting='none';

if ~isfield(Anat,'ctx'),Anat.ctx='std';end
switch Anat.ctx
    case 'vinf'
Rkey=10;
Xkey=-115;
Ykey0=100;
    case 'inf'
Rkey=8;
Xkey=-105;
Ykey0=100;
    case 'std'
Rkey=5;
Xkey=-75;
Ykey0=60;
end
 if ~isfield(params,'radius'), params.radius=repmat(5,[Nroi,1]);end
 radius= params.radius;
 ROI.radius = repmat(radius,[Nroi,1]);
if ~isfield(ROI,'radius'), ROI.radius=repmat(Rkey,[Nroi,1]);end


if ~isfield(Anat,'alpha'),Anat.alpha=1;end


%% Display figures 
f=figure('Position',[100,100,1550,750]);

ROI2.radius = repmat(radius,[Nroi,1]);
ROI2.color=Cmap(Modules,:);
ROI2.coord=ROIxyz;
if isfield(ROI2,'Network'),ROI2=rmfield(ROI2,'Network');end
ROI2.Network(:,1)=1:size(ROI2.coord,1);
ROI2.Network(:,2)=ones(size(ROI2.coord,1),1);

subplot(2,3,[2:3],'Position',[.45,0.455,.53,.48])
Anat.view='lat';Anat.alpha=1;Draw_ROIs_on_Cortex(Anat,ROI2);

subplot(2,3,[5:6],'Position',[.45,0.005,.53,.48])
Anat.view='med';Anat.alpha=1;Draw_ROIs_on_Cortex(Anat,ROI2);

subplot(2,3,[1,4],'Position',[.075,0.025,.35,.9])
Anat.view='dorsal';Draw_ROIs_on_Cortex(Anat,ROI2);
light('Position',[0,100,100],'Style','local');
set(gcf,'Color','w')
title(strrep(IM.name,'_',' '))

% axis on
% xlabel('x');ylabel('y');zlabel('z')

%% Add key ROIs
for j=1:length(Cmap)
    foci.location(j,:)=[Xkey,Ykey0-(j-1)*13,30];
    foci.radius(j)=Rkey;
    foci.color(j,:)=Cmap(j,:);
    text(Xkey-7,Ykey0-(j-1)*13,30,IM.Nets(j,1),...
        'HorizontalAlignment','right','FontName','Arial','FontSize',15);
end
foci.lighting=Params.lighting;
Draw_Foci(foci,ROI.Nfaces);
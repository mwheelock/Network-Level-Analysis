function foci=AdjustFoci(foci,meshL,meshR,side,ctx)
%
% This function adjusts Foci locations to match the surface of the mesh
% in case it is inflated or rotated or otherwise adjusted.
%
% New algorithm 8/16/17, to be fixed....
% 1. deal with std vs inf or vinf first
% 2. determine home hemisphere for ROIs
% 3. keep track of full transformation on mesh surface for ROIs
% 4. apply to ROIs L and R sep given what happened to hem



%% First, adjust to match the ctx being imaged
if ~exist('params','var'),params=struct;end

%% Determine home hemisphere
pos=foci.location;

[IDXl,dl]=knnsearch(meshL.nodes,pos);
[IDXr,dr]=knnsearch(meshR.nodes,pos);

focL=dl<dr;
focR=dr<=dl;


%% Choose inflation and adjust foci to match ctx
switch ctx
    case 'std'
        Lnodes=meshL.nodes;
        Rnodes=meshR.nodes;
    case 'inf'
        Lnodes=meshL.Inodes;
        Rnodes=meshR.Inodes;
    case 'vinf'
        Lnodes=meshL.VInodes;
        Rnodes=meshR.VInodes;
end
dxL=Lnodes-meshL.nodes;
dxR=Rnodes-meshR.nodes;

pos(focL,:)=pos(focL,:)+dxL(IDXl(focL),:);
pos(focR,:)=pos(focR,:)+dxR(IDXr(focR),:);


%% Small adjustment to separate hemispheres
pos(focL,1)=pos(focL,1)-max(Lnodes(:,1));
pos(focR,1)=pos(focR,1)-min(Rnodes(:,1));

Lnodes(:,1)=Lnodes(:,1)-max(Lnodes(:,1));
Rnodes(:,1)=Rnodes(:,1)-min(Rnodes(:,1));


%% Rotate if necessary
if (strcmp(side,'lat') || strcmp(side,'med'))
dy=-5;
    % rotate right hemi around and move to position for visualization
    cmL=mean(Lnodes,1);
    cmR=mean(Rnodes,1);
    rm=RotationMatrix('z',pi);
    
    % Rotate
    switch side
        case 'lat'
            Rnodes=(Rnodes-(repmat(cmR,size(Rnodes,1),1)))*rm +...
                (repmat(cmR,size(Rnodes,1),1));             
            pos(focR,:)=(pos(focR,:)-(repmat(cmR,sum(focR),1)))*rm +...
                (repmat(cmR,sum(focR),1));             
        case 'med'
            Lnodes=(Lnodes-(repmat(cmL,size(Lnodes,1),1)))*rm +...
                (repmat(cmL,size(Lnodes,1),1));       
            pos(focL,:)=(pos(focL,:)-(repmat(cmL,sum(focL),1)))*rm +...
                (repmat(cmL,sum(focL),1));              
    end
    pos(focR,1)=pos(focR,1)+(cmL(:,1)-cmR(:,1));    % Shift over to same YZ plane
    pos(focR,2)=pos(focR,2)-max(Rnodes(:,2))+min(Lnodes(:,2))+dy;
end


%% Put back in output
foci.location=pos;
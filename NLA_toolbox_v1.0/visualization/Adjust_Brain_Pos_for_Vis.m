function [Rnodes,Lnodes]=Adjust_Brain_Pos_for_Vis(meshL,meshR,params,side,ctx)
%
% This function takes as inputs the meshes for the 2 hemispheres, along
% with the inflation and perspective information for plotting.
% The function positions the meshes for ideal visualizations.

%% Parameters and Initialization
if ~isfield(params,'info');params.info=struct;params.info.acq='t';end


%% Choose inflation
switch ctx
    case 'std'
        Lnodes=meshL.nodes;
        Rnodes=meshR.nodes;
        
        % If volume is not in transverse orientation, put it there
        switch params.info.acq(1)
            case 's'
                Nln=size(Lnodes,1);
                temp=cat(1,Lnodes,Rnodes);
                temp=RotateCap(temp,[-90,0,90]);
                Lnodes=temp(1:Nln,:);
                Rnodes=temp((Nln+1):end,:);
            case 'c'
                % not yet supported
        end

    case 'inf'
        Lnodes=meshL.Inodes;
        Rnodes=meshR.Inodes;
    case 'vinf'
        Lnodes=meshL.VInodes;
        Rnodes=meshR.VInodes;
end



%% Small adjustment to separate hemispheres
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
        case 'med'
            Lnodes=(Lnodes-(repmat(cmL,size(Lnodes,1),1)))*rm +...
                (repmat(cmL,size(Lnodes,1),1));              
    end
    Rnodes(:,1)=Rnodes(:,1)+(cmL(:,1)-cmR(:,1));    % Shift over to same YZ plane
    Rnodes(:,2)=Rnodes(:,2)-max(Rnodes(:,2))+min(Lnodes(:,2))+dy;
end
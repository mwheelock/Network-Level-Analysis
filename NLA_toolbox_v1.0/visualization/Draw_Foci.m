function Draw_Foci(foci,faces)

% This function draw spheres of a specified color and radius at the
% location given by the structure foci. Fields include: location, color,
% radius. It is assumed that coordinate space and voxel-index space are
% exactly aligned. foci.color are rows of RGB values [0,1].

if ~exist('faces','var'), faces=100;end
if ~isfield(foci,'lighting'),foci.lighting='phong';end

[x,y,z]=sphere(faces);
for j=1:size(foci.location,1)
    
    rad=foci.radius(j);
    
    hh=patch(surf2patch(rad*x+foci.location(j,1),...
        rad*y+foci.location(j,2),rad*z+foci.location(j,3)),...
        'EdgeColor',foci.color(j,:),'FaceColor',foci.color(j,:),...
        'EdgeAlpha',0);
    
    %    set(hh,'FaceLighting','phong','AmbientStrength',0.02);
    %   set(hh,'SpecularStrength',0.01,'AmbientStrength',.1);
    %    set(hh,'FaceLighting','flat','AmbientStrength',0.02);
    %    set(hh,'FaceLighting','gouraud','AmbientStrength',0.02);
    %    set(hh,'FaceLighting','none','AmbientStrength',0.02);
    set(hh,'FaceLighting',foci.lighting,'AmbientStrength',0.02);
    
end

axis image
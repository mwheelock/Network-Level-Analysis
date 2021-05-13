function Draw_ROIs_Through_Cortex_3_Views(Anat,roi,Conn,params)

% This function displays hemispheres with ROI spheres and sticks connecting
% them in 3 views.

Px=0;
Py=-110;
Pz=-125;
Txt='+';
Ctxt=[1,0,0];


disp('<<< Drawing ROIs and sticks on transparent brains, please be patient')
figure('Color','w','Position',[50 300 1650 750])
a1=subplot(1,3,1,'Position',[0.14,0.1,0.24,0.84]);  % Posterior View [left,bottom,w,h]
Draw_ROIs_Through_Cortex(Anat,roi,Conn,params);
view([0,0])
%pause(1); 
a2=subplot(1,3,2,'Position',[0.4,0.1,0.22,0.82]);  % Dorsal View
%pause(1);
a2=subplot(1,3,2,'Position',[0.4,0.1,0.22,0.82]);  % Dorsal View
Draw_ROIs_Through_Cortex(Anat,roi,Conn,params);
view([0,90])
% title(params.group{g})
a3=subplot(1,3,3,'Position',[0.64,0.08,0.27,0.87]);  % Left Lateral View
Draw_ROIs_Through_Cortex(Anat,roi,Conn,params);
view([-90,0])
text(Px,Py,Pz,Txt,'Color',Ctxt)

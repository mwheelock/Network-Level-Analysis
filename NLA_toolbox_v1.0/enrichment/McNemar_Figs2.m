function EWsig=McNemar_Figs2(dataOut,IM,B,params)

% Plot network pairs meeting significance

if ~exist('B','var'),B=1;end
Pth=0.05./B;
PthMcN=0.05;
[NNidx,Tidx,TNidx]=IM2idx(IM);
Npairs=length(TNidx);
Nnets=max(IM.key(:,2));
NetKey=[[1:Nnets]',[1:Nnets]'];
%cMap=[1,1,1;0,0,0];

%% Significant stats by network block
gt1 = char(strcat({'Enrichment '},params.group{1})); % name of enriched group1/time1
gt1M = char(strcat({'Enriched & McNemar '}, params.group{1}));
gt2 = char(strcat({'Enrichment '}, params.group{2}));
gt2M = char(strcat({'Enriched & McNemar '}, params.group{2}));
g1g2 = char(strcat({'Enriched '},params.group{1},{' & '},params.group{2})); % both groups enriched
unicorn = char(strcat({'Enriched & McNemar '}, params.group{1},{' & '},params.group{2})); % Enrichment and McNemar
cMap2=[1,1,1;0.5,0.75,1;0,0,1;...% none, EN12(noMnC),EN12(McN)
    1,0.6,0.6;1,0,0;... %EN24(noMcN),EN24(McN)
    0,1,0;1,0,1]; %EN12,24(noMcN),EN12,24(McN)
Xpos=8;
FS=16;
Tcol='k';
Fcol='w';
E12=(squeeze(dataOut.Chi_EWpval(:,:,1))<=Pth).*...
    (squeeze(dataOut.HGppEW(:,:,1))<=Pth);
E24=(squeeze(dataOut.Chi_EWpval(:,:,2))<=Pth).*...
    (squeeze(dataOut.HGppEW(:,:,2))<=Pth);
McN=dataOut.McNemar_PPEW<=Pth;
a=(E12.*(-E24+1).*(-McN+1)).*1+...% E12, no E24, no McN
    (E12.*(-E24+1).*(McN)).*2+...% E12, no E24,  McN;
    (E24.*(-E12+1).*(-McN+1)).*3+...% E24, no E12, no McN;
    (E24.*(-E12+1).*(McN)).*4+...% E24, no E12,  McN;
    (E12.*(E24).*(-McN+1)).*5+...% E12, E24, no McN;
    (E12.*(E24).*(McN)).*6; % E12, E24, McN;
figure('Color',Fcol,'Position',[100,100,900,600]);
Matrix_Org3(a,NetKey,0.5,[0,6],IM.cMap,0,cMap2);
hold on;
rectangle('Position',[Xpos,0.5,0.5,0.5],'FaceColor',cMap2(2,:))
text(Xpos+1,0.75,gt1,'Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,1.5,0.5,0.5],'FaceColor',cMap2(3,:))
text(Xpos+1,1.75,gt1M,'Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,2.5,0.5,0.5],'FaceColor',cMap2(4,:))
text(Xpos+1,2.75,gt2,'Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,3.5,0.5,0.5],'FaceColor',cMap2(5,:))
text(Xpos+1,3.75,gt2M,'Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,4.5,0.5,0.5],'FaceColor',cMap2(6,:))
text(Xpos+1,4.75,g1g2,'Color',Tcol','FontSize',FS)
rectangle('Position',[Xpos,5.5,0.5,0.5],'FaceColor',cMap2(7,:))
text(Xpos+1,5.75,unicorn,'Color','k','FontSize',FS)
% set(gcf,'Color','k')

EWsig = a;
% 1) group 1 enriched 
% 2) group 2 enriched
% 3) group 1 mcnemar
% 4) group 2 mcnemar
% 5) shared enrichment
% 6) enrichment and mcnemar


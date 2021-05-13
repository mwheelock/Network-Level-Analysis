function EWsig=Enrichment_Figs2(dataOut,IM,B)

% Plot network pairs meeting significance

if ~exist('B','var'),B=1;end
Pth=0.05./B;
PthMcN=0.05;
[NNidx,Tidx,TNidx]=IM2idx(IM);
Npairs=length(TNidx);
Nnets=max(IM.key(:,2));
NetKey=[[1:Nnets]',[1:Nnets]'];
cMap=[1,1,1;0,0,0];

%% Significant stats by network block
Xpos=8;
Tcol='k';
Fcol='w';
% 
% cMap2=[1,1,1;0.5,0.75,1;0,0,1;...% none, EN12(noMnC),EN12(McN)
%     1,0.5,0;1,0,0;... %EN24(noMcN),EN24(McN)
%     0,1,0;1,0,1]; %EN12,24(noMcN),EN12,24(McN)
cMap2=[1,1,1;0,0,0];
dataOut.Enrich_sig=(squeeze(dataOut.Chi_EWpval(:,:,1))<=Pth).*...
    (squeeze(dataOut.HGppEW(:,:,1))<=Pth);

figure('Color',Fcol,'Position',[100,100,900,600]);
Matrix_Org3(dataOut.Enrich_sig,NetKey,0.5,[0,1],IM.cMap,0,cMap2);


EWsig=dataOut.Enrich_sig;


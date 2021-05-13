% Script for visualizing Enrichment results

%% Visualize thresholded network matrix

% load enrichment results and corresponding IM structure
load('name of structure containing dataOut');
load('name of IM structure');
fc = load('name of fc data'); % rvals, not zvals
%fc=fc.rmat;
fc=FisherR2Z(fc.rmat);   
fc=FisherZ2R(fc);

% visualize matrix of network pairs meeting signficance threshold
B = 1; % set your bonferroni threshold (.05/B)
Perm_1tp_Figs2(dataOut,IM,B,params);
EWsig=Enrichment_Figs2(dataOut,IM,B);

%% Visualize significant networks on the brain

%load('Conte69_on_TT_32k.mat','Anat'); %(Talairach)
load('MNI_coord_meshes_32k','MNIl','MNIr'); %(MNI 152 in MNI)
%load('CH2_MNI_withCB','MNIl','MNIr'); %(MNI Colin Holmes with cerebellum)
Anat.CtxL=MNIl;Anat.CtxR=MNIr;
Anat.alpha = 0.75 ;
Anat.ctx='std'; % typically 'inf' for adult mesh
clear MNIl MNIr

% 'Balls on Brains' - coloring by fcBx correlation only
params.ScaleRadius=0; % 1 = set ROI size to reflect degree 
params.roiradius=5; % 5 is recommended for adult atlas space

% Matrix by behavior % 
[r,c]=ind2sub(size(EWsig),find(EWsig(:)));
RC=[r,c];
g =1;

% RC= [7,3;10,2;10,8;11,3;13,8; 14,3; 15,8; 15, 10; 16,8];

% Visualize each brain network pair separately 
% NOTE: press space bar to continue
for j=1:size(RC,1)
pairs{j}=View_Clear_Brain_ROI_Corr(Anat,IM,RC(j,:),[],...[Row,Col]
    squeeze(dataOut.pval(:,:,1)),squeeze(dataOut.rho(:,:,1)),...
    0.05,params,g);
% pause;
%close; 

mask=zeros(size(IM.key,1));
mask(IM.key(:,2)==RC(j,1),IM.key(:,2)==RC(j,2))=1;
pj=dataOut.pval(mask==1);pj(pj==0)=1;
rj=dataOut.rho(mask==1);

figure('Color','w','Units','Normalized','Position',[0.35,0.56,0.47,0.35]);
histogram(rj(rj~=0),[-1:1e-2:1],'FaceColor','k');hold on
histogram(rj(((pj<0.05).*(rj<0))>0),[-1:1e-2:1],'FaceColor','b');
histogram(rj(((pj<0.05).*(rj>0))>0),[-1:1e-2:1],'FaceColor','r');
set(gca,'YScale','log');yl=ylim;ylim([0.5,yl(2)])
xlabel([params.type,' Correlation']);ylabel('Number ROI pairs')
text(0.5,5,[{['Total hits: ',num2str(100*sum(pj<0.05)/length(pj)),'%']};...
    {['p(\chi^2)=',num2str(dataOut.Chi_EWpval(RC(j,1),RC(j,2)))]};...
    {['p(HGp)=',num2str(dataOut.HGppEW(RC(j,1),RC(j,2)))]}])
title([IM.Nets{RC(j,1)},'-',IM.Nets{RC(j,2)},...
    ' Univariate fc-behavior correlations']) 
end
 
%Visualize all network pairs on one brain (not recommended)
View_Clear_Brain_ROI_Corr_Multi(Anat,IM,RC,[],...[Row,Col]
    squeeze(dataOut.pval(:,:,1)),squeeze(dataOut.rho(:,:,1)),...
    0.05,params,g);

%Visualize each network pair with direction of fc
for j=1:length(RC)
    pairsOut{j}=View_Clear_Brain_ROI_Corr_fc(Anat,IM,RC(j,:),...
        dataOut,fc1,0.05,params);
end


%% Looking more deeply at specific ROI-pairs within implicated network pairs

% The pairs cell arrays above list the specific ROI pairs 
% that drive the enrichment results.
% Use them to further explore the specific brain-behavior (fcBx)
% relationships.

% Scatter plot for a specific ROI pair within a specific network pair
fnPair=1; % functional network pair within pairsout array
pairNum=20;% ROI pair. these indices reference the entire fc matrix.

ROI1=pairs{1,fnPair}{1,1}(pairNum,1);
ROI2=pairs{1,fnPair}{1,1}(pairNum,2);
rp=squeeze(fc(ROI1,ROI2,:));
% yaxismin = min(Bx)-std(Bx);
% yaxismax = max(Bx)+std(Bx);

figure('Color','w');
plot(rp,Bx,'ok');xlabel('fc [z(r)]');ylabel([params.BxName]);
xlim([-1 1]);
switch params.type
    case 'Pearson'
lsline
title([{[IM.Nets{RC(fnPair,1)},'-',IM.Nets{RC(fnPair,2)}]};...
    {['ROI ',num2str(ROI1),' - ROI ',num2str(ROI2)]};...
    {[num2str(IM.ROIxyz(ROI1,:)),';',num2str(IM.ROIxyz(ROI2,:))]};...
    {['r=',num2str(dataOut.rho(ROI1,ROI2)),...
    '; p=',num2str(dataOut.pval(ROI1,ROI2))]}])

    case 'Spearman'
title([{[IM.Nets{RC(fnPair,1)},'-',IM.Nets{RC(fnPair,2)}]};...
    {['ROI ',num2str(ROI1),' - ROI ',num2str(ROI2)]};...
    {[num2str(IM.ROIxyz(ROI1,:)),';',num2str(IM.ROIxyz(ROI2,:))]};...
    {['rho=',num2str(dataOut.rho(ROI1,ROI2)),...
    '; p=',num2str(dataOut.pval(ROI1,ROI2))]}])
end

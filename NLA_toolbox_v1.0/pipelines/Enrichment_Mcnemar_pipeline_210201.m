%% Enrichment and Mcnemar
% This code runs enrichment on 2 independent groups. 
% McNemar chi-square is used to compare network pairs between groups. 

%% Load data

% load IM structure
load('filepath to IM structure');

% load connectivity data
rmat1=load('filepath to fc dataset 1');
rmat2=load('filepath to fc dataset 12');
% re-order rmats to the network order in IM (if already ordered - skip)
rmat1ordered= rmat1(IM.order,IM.order,:);
rmat2ordered= rmat2(IM.order,IM.order,:);

% load behavior data
behavdata=importdata('filepath to behavioral data');

%% plot average fc organized by networks

% if r-values
    rmatAve=mean(FisherR2Z(rmat),3); 
% if r-values
    rmatAve=mean(rmat,3);

% % load anat (surface file) % %
load('MNI152nl_on_TT_coord_meshes_32k','MNIl','MNIr'); % MNI152 template in TT
%load('MNI_coord_meshes_32k','MNIl','MNIr'); % MNI152 in MNI
Anat.CtxL=MNIl;Anat.CtxR=MNIr;
Anat.alpha = 0.25; % how see through is mesh?
clear MNIl MNIr

% % Visualize rmat organized by IM file % %
Matrix_Org3(rmatAve(IM.order,IM.order),...
    IM.key,10,[-0.3,0.3],IM.cMap,1);

% % Visualize IM Model on Brain with Network Names and Colors 
params.radius = 4; %specify ROI size for display on brain surface
Anat.ctx='std'; %typically use 'inf' for adult mesh, mesh with CB is 'std'
View_ROI_Modules(IM,Anat,IM.ROIxyz,params);

% % plot ROI count in each network % %
figure;histogram(IM.key(:,2));title(strrep(IM.name,'_',' '))
set(gca,'XTick',[1:max(IM.key(:,2))],'XTickLabel',IM.Nets);ylabel('Nrois')

% specify names of groups for McNemar data structure
group1 = 'group1';
group2 = 'group2';
Behavename = 'behavior name';

% spedify behaviors to test in the model
behavior1 = Bx1;
behavior2= Bx2;

% make cell arrays to input into McNemar
Bx = {behavior1 behavior2};
fc = {rmat1ordered rmat2ordered};
%% 4. Enrichment & McNemar Analyses

params.np = 10000; % number permutations
params.type='Pearson'; % change to 'Spearman' if Bx does not meet normality
params.tp = 1; % params.tp=2 goes into a time point cohort select code specific to IBIS data
params.Pmax=0.05; % sets the ROI-level p value threshold
params.nnPmax=0.05; % 
params.noUS=1;
params.B = 1; %1 = .05, 2=0.025, and so forth, recommend to run at low threshold
params.group = {group1 group2};
params.BxName = Behavename;
dataIN.fc=fc; 
dataIN.Bx=Bx;
dataIN.IM=IM;
dataIN.params=params;
params.fn=strcat(Behavename,num2str(params.np),'_',num2str(params.Pmax),'_',IM.name,'.mat'); 
            
% Run Enrichment and McNemar with a p threshold
[dataOut]=fcBx_Chi_McNemar_TwoGroups(dataIN.fc,dataIN.Bx,dataIN.IM,dataIN.params);

% save data
cd('filepath to output directory');
save(params.fn,'dataOut','dataIN','params');
clear dataOut dataIN

% save data
cd('filepath to output directory');
save(params.fn,'dataOut','dataIN','params');
clear dataOut dataIN

%% Extract Chi-square and McNemar test p-values
RC = [5,4;10,2;8,2;12,5]; % specify row,column of network statistic to extract

Nnets=length(IM.Nets);
sizenets=dataOut.Np;
for j=1:length(dataOut.TNidx);
    nsize(dataOut.TNidx(j))=sizenets(j);
end
netsize=reshape(nsize,Nnets,Nnets,[]);

params.roiradius=5;
params.toDraw=0;
for j=1:length(RC)    
pairsOut{j}=View_Clear_Brain_ROI_Corr_fc2(Anat,dataIN.IM,RC(j,:),...
    dataOut,dataIN.fc,0.05,params);
end
for j=1:size(RC,1)
    for g=1:length(params.group);
        name=params.group(g);name=name{1,1};
NetStats.nets{j} = strcat(dataIN.IM.Nets{RC(j,1)},'-',dataIN.IM.Nets{RC(j,2)});
NetStats.(name).Chipvals(j) =dataOut.Chi_EWpval(RC(j,1),RC(j,2),g);
NetStats.(name).Chistats(j) =dataOut.Chi_stats(RC(j,1),RC(j,2),g); 
NetStats.McNemarpvals(j) = dataOut.McNemar_PPEW(RC(j,1),RC(j,2)); 
NetStats.(name).Nethits(j)=size(pairsOut{1,j}{1,1},g);
NetStats.(name).N(j) = netsize(RC(j,1),RC(j,2));
NetStats.(name).phieffectsize(j)=sqrt(dataOut.Chi_stats(RC(j,1),RC(j,2),g)/(netsize(RC(j,1),RC(j,2))));
NetStats.McNemar_stats(j) = dataOut.McNemar_stats(RC(j,1),RC(j,2)); 
NetStats.McNemarphieffectsize(j)=sqrt(dataOut.McNemar_stats(RC(j,1),RC(j,2))/(netsize(RC(j,1),RC(j,2))));
    end
end

%% Enrichment and McNemar Visualizations

% replot data at more stringent threshold if desired
% BrBx correlations will still be p<0.05, but group level significance can be bonferroni (B) corrected if desired
B=1; % 1=p<0.05; 2=p<0.025; etc.

%re-plot all figures if desired%
Perm_Mixed_LCS_Figs2(dataOut,IM,B,params);

% re-plot the triangle results only and generate vector for making brain visualizaitons %
EWsig=McNemar_Figs2(dataOut,IM,B,params); %experiment wide significance

%% make behavior histogram

% Example 1 %
d1=dataIN.Bx{1,1};
[test1 test2] = hist(d1);
d2=dataIN.Bx{1,2};
mid1=min(d1); mid2=min(d2); mad1=max(d1); mad2=max(d2);
if mid1>mid2 ;usemin = mid2; else usemin = mid1; end
if mad1>mad2; usemax = mad1; else usemax = mad2; end
myBins = linspace(usemin,usemax,15); % specify your desired number of bins
[y1]=hist(d1,myBins); 
[y2]=hist(d2,myBins);
figure(1)
set(gcf,'color','w');
bar(myBins, [y1;y2]');
% change your color, width, font size manually using edit>figure properties

% Example 2 %
d1=dataIN.Bx{1,1};
d2=dataIN.Bx{1,2};
figure 
set(gcf,'color','w');
histogram(d1,usemin:2:usemax,'facecolor', 'red','facealpha',.4, 'edgecolor','none')
hold on
histogram(d2,usemin:2:usemax,'facecolor', 'blue','facealpha',.4, 'edgecolor','none')
box off
axis tight
legend('VPT','Term','location','northeast');

%% Plot data on the brain
RC = [5,4;7,5;15,5]; % specify row,column of network statistic to extract
params.toDraw=1;
for g=1:length(params.group);
View_Clear_Brain_ROI_Corr_Multi(Anat,dataIN.IM,RC,[],...[Row,Col]
    squeeze(dataOut.pval(:,:,g)),squeeze(dataOut.rho(:,:,g)),...
    0.05,params,g);
end

for j=1:length(RC)    
pairsOut{j}=View_Clear_Brain_ROI_Corr_fc2(Anat,dataIN.IM,RC(j,:),...
    dataOut,dataIN.fc,0.05,params);
end

%% Scatter plot for a specific ROI pair within a specific network pair

g = 1; % specify which group to plot (1 or 2)
fnPair=1; % functional network pair within pairsOut cell array
pairNum=22;% ROI network pair within pairsOut

ROI1=pairsOut{1,fnPair}{1,g}(pairNum,1);
ROI2=pairsOut{1,fnPair}{1,g}(pairNum,2);
rp=squeeze(fc{g}(ROI1,ROI2,:));


figure('Color','w');
plot(rp,Bx{g},'ok');xlabel('fc [z(r)]');ylabel([params.BxName]);
xlim([-1 1]);
switch params.type
    case 'Pearson'
lsline
title([{[IM.Nets{RC(fnPair,1)},'-',IM.Nets{RC(fnPair,2)}]};...
    {['ROI ',num2str(ROI1),' - ROI ',num2str(ROI2)]};...
    {[num2str(IM.ROIxyz(ROI1,:)),';',num2str(IM.ROIxyz(ROI2,:))]};...
    {['r=',num2str(dataOut.rho(ROI1,ROI2,g)),...
    '; p=',num2str(dataOut.pval(ROI1,ROI2,g))]}])

    case 'Spearman'
        lsline
title([{[IM.Nets{RC(fnPair,1)},'-',IM.Nets{RC(fnPair,2)}]};...
    {['ROI ',num2str(ROI1),' - ROI ',num2str(ROI2)]};...
    {[num2str(IM.ROIxyz(ROI1,:)),';',num2str(IM.ROIxyz(ROI2,:))]};...
    {['rho=',num2str(dataOut.rho(ROI1,ROI2,g)),...
    '; p=',num2str(dataOut.pval(ROI1,ROI2,g))]}])
end

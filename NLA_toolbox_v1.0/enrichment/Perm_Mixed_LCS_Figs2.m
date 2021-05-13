function Perm_Mixed_LCS_Figs2(dataOut,IM,B,params)

% This function is the same as Perm_Mixed_LCS_Figs except that it does not
% plot the block-wise FPR values. Also, it plots the FPR-p-vals for both
% Chi-squared and HG to see dynamic range
if ~exist('B','var'),B=1;end
Pth=0.05./B;
PthMcN=0.05;
[NNidx,Tidx,TNidx]=IM2idx(IM);
Npairs=length(TNidx);


Nnets=max(IM.key(:,2));
NetKey=[[1:Nnets]',[1:Nnets]'];
cMap=[1,1,1;0,0,0];
cMap2=[1,1,1;1,0,0;0,0,1;1,0,1];

%% McNemar
figure('Color','w','Position',[100,100,1000,800]);
subplot(2,2,1)
plot([0,1],-log10([dataOut.McN_EWth,dataOut.McN_EWth]),'r','LineWidth',2);
hold on;xlabel('b/(b+c)');ylabel('-log_1_0(p)');
title('McNemar-asymptotic');
L10p=-log10(dataOut.McNemar_P(TNidx));
plot(dataOut.McNemar_b(TNidx)./...
    (dataOut.McNemar_b(TNidx)+dataOut.McNemar_c(TNidx)),L10p,'k*');
if max(L10p)>8, axis([0,1,0,ceil(max(L10p))]),else axis([0,1,0,8]);end
subplot(2,2,3);
plot([0,1],-log10([PthMcN,PthMcN]),'b','LineWidth',2);hold on
L10pp=-log10(dataOut.McNemar_PPEW(TNidx));
plot(dataOut.McNemar_b(TNidx)./...
    (dataOut.McNemar_b(TNidx)+dataOut.McNemar_c(TNidx)),L10pp,'k*');
xlabel('b/(b+c)');ylabel('-log_1_0(p)');title('McNemar - permuted');
if max(L10pp)>3, axis([0,1,0,ceil(max(L10p))]),else axis([0,1,0,3]);end
subplot(2,2,2)
loglog(dataOut.Emp_FDR_McN(:,1),dataOut.Emp_FDR_McN(:,2),'k');hold on;
loglog([1e-20,1],[PthMcN,PthMcN],'b','LineWidth',2)
loglog([dataOut.McN_EWth,dataOut.McN_EWth],[1e-5,1],'r','LineWidth',2)
loglog(dataOut.McNemar_P(TNidx),dataOut.McNemar_PPEW(TNidx),'*k');
title('McNemar p-values');xlabel('Asymptotic');
ylabel('Permutation-based FPR')
loglog([0.05,0.05],[1e-5,1e0],'k')
loglog([0.05,0.05]./Npairs,[1e-5,1e0],'g')
axis([10^(-ceil(-log10(min(dataOut.McNemar_P(TNidx))))-2),1,...
    10^(-ceil(-log10(min(dataOut.McNemar_PPEW(TNidx))))-1),1])
subplot(2,2,4);
Matrix_Org3(-log10(dataOut.McNemar_PPEW),NetKey,0.5,[0,2],IM.cMap,0,parula);
title('McNemar')
[r,c]=ind2sub(size(dataOut.McNemar_PPEW),find(dataOut.McNemar_PPEW<Pth));
hold on
plot(c,r,'*k');colorbar


%% Enrichment 
a=reshape(dataOut.Chi_pval0,[],2);b=reshape(dataOut.Chi_EWpval,[],2);
c=reshape(dataOut.HGp,[],2);d=reshape(dataOut.HGppEW,[],2);
for j=1:2
    if j==1
        tempname = params.group{1};
    else
        tempname = params.group{2};
    end
figure('Color','w','Position',[100,100,1400,800]);
subplot(2,3,1)
Matrix_Org3(squeeze(dataOut.rho(:,:,j)),IM.key,10,[-0.3,0.3],IM.cMap,0,jet);
title([tempname,' ',params.type])
subplot(2,3,4)
Matrix_Org3(squeeze(dataOut.thresholdedROIpairs(:,:,j)),IM.key,10,...
    [0,1],IM.cMap,0,cMap);
title([tempname,' ',params.type,' p<0.05'])

subplot(2,3,2)
loglog(squeeze(dataOut.Emp_FDR_CS(:,1,j)),...
    squeeze(dataOut.Emp_FDR_CS(:,2,j)),'k');hold on;
loglog([1e-20,1],[Pth,Pth],'b','LineWidth',2)
loglog([squeeze(dataOut.Chi_EWth(j)),squeeze(dataOut.Chi_EWth(j))],...
    [1e-5,1],'r','LineWidth',2)
loglog(squeeze(a(TNidx,j)),squeeze(b(TNidx,j)),'*k');
title('\chi^2 p-values');xlabel('Asymptotic');
ylabel('Permutation-based FPR')
loglog([0.05,0.05],[1e-5,1e0],'k')
loglog([0.05,0.05]./Npairs,[1e-5,1e0],'g')
axis([10^(-ceil(-log10(min(a(TNidx,j))))-2),1,...
    10^(-ceil(-log10(min(b(TNidx,j))))-1),1])

subplot(2,3,5);
foo=-log10(dataOut.Chi_EWpval(:,:,j));
Matrix_Org3(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx)),2],IM.cMap,0,parula(1000));
[r,col]=ind2sub(size(squeeze(dataOut.Chi_EWpval(:,:,j))),...
    find(squeeze(dataOut.Chi_EWpval(:,:,j))<Pth));
hold on
plot(col,r,'*k');%colorbar


subplot(2,3,3)
loglog(squeeze(dataOut.Emp_FDR_HG(:,1,j)),...
    squeeze(dataOut.Emp_FDR_HG(:,2,j)),'k');hold on;
loglog([1e-20,1],[Pth,Pth],'b','LineWidth',2)
loglog([squeeze(dataOut.HG_EWth(j)),squeeze(dataOut.HG_EWth(j))],...
    [1e-5,1],'r','LineWidth',2)
loglog(squeeze(c(TNidx,j)),squeeze(d(TNidx,j)),'*k');
title('Hypergeometric p-values');xlabel('Asymptotic');
ylabel('Permutation-based FPR')
loglog([0.05,0.05],[1e-5,1e0],'k')
loglog([0.05,0.05]./Npairs,[1e-5,1e0],'g')
axis([10^(-ceil(-log10(min(c(TNidx,j))))-2),1,...
    10^(-ceil(-log10(min(d(TNidx,j))))-1),1])

subplot(2,3,6)
foo=-log10(dataOut.HGppEW(:,:,j));
Matrix_Org3(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx)),2],IM.cMap,0,parula(1000));
[r,col]=ind2sub(size(squeeze(dataOut.HGppEW(:,:,j))),...
    find(squeeze(dataOut.HGppEW(:,:,j))<Pth));
hold on
plot(col,r,'*k');%colorbar
end

% Plot within group enrichment in addition to McNemar?
if ~isfield(params,'PlotEnrich'), params.PlotEnrich='WithEnrich';end
PlotEnrich=params.PlotEnrich;
if ~exist('PlotEnrich','var'), PlotEnrich=WithEnrich; end
switch PlotEnrich
    case 'NoEnrich'
            cMap2=[1,1,1;01,1,1;0,0,1;...% none, EN12(noMnC),EN12(McN)
            1,1,1;1,0,0;... %EN24(noMcN),EN24(McN)
            0,1,0;1,0,1]; %EN12,24(noMcN),EN12,24(McN)
    case 'WithEnrich'
            cMap2=[1,1,1;0.5,0.75,1;0,0,1;...% none, EN12(noMnC),EN12(McN)
            1,0.5,0;1,0,0;... %EN24(noMcN),EN24(McN)
            0,1,0;1,0,1]; %EN12,24(noMcN),EN12,24(McN)
end
        
%% Significant stats by network block
gt1 = char(strcat({'Enrichment '},params.group{1})); % name of enriched group1/time1
gt1M = char(strcat({'Enriched & McNemar '}, params.group{1}));
gt2 = char(strcat({'Enrichment '}, params.group{2}));
gt2M = char(strcat({'Enriched & McNemar '}, params.group{2}));
g1g2 = char(strcat({'Enriched '},params.group{1},{' & '},params.group{2})); % both groups enriched
unicorn = char(strcat({'Enriched & McNemar '}, params.group{1},{' & '},params.group{2})); % Enrichment and McNemar
% cMap2=[1,1,1;0.5,0.75,1;0,0,1;...% none, EN12(noMnC),EN12(McN)
%     1,0.5,0;1,0,0;... %EN24(noMcN),EN24(McN)
%     0,1,0;1,0,1]; %EN12,24(noMcN),EN12,24(McN)
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

%% Double check Network-pair size is not a confound
% a=reshape(dataOut.Chi_EWpval,[],2);
% b=reshape(dataOut.HGppEW,[],2);
a=reshape(dataOut.Chi_pval0,[],2);
b=reshape(dataOut.HGp,[],2);
figure;
subplot(2,2,1)
plot(dataOut.Np,-log10(a(TNidx,1)),'ko');hold on
plot(dataOut.Np,-log10(b(TNidx,1)),'r+');hold on
xlabel('Number of ROI pairs within network pair')
ylabel('Asymptotic p-value')
title(params.group{1})
[r,p]=corr(dataOut.Np,a(TNidx,1));
text(700,7,['CS r = ',num2str(r)])
text(700,6,['CS p = ',num2str(p)])
[r,p]=corr(dataOut.Np,b(TNidx,1));
text(700,5,['HG r = ',num2str(r)])
text(700,4,['HG p = ',num2str(p)])
% axis([0,1000,1e-4,1])

subplot(2,2,2)
plot(dataOut.Np,-log10(a(TNidx,2)),'ko');hold on
plot(dataOut.Np,-log10(b(TNidx,2)),'r+');hold on
xlabel('Number of ROI pairs within network pair')
ylabel('Asymptotic p-value')
title(params.group{2})
[r,p]=corr(dataOut.Np,a(TNidx,2));
text(700,5,['CS r = ',num2str(r)])
text(700,4,['CS p = ',num2str(p)])
[r,p]=corr(dataOut.Np,b(TNidx,2));
text(700,3,['HG r = ',num2str(r)])
text(700,2,['HG p = ',num2str(p)])
% axis([0,1000,1e-4,1])

a=reshape(dataOut.Chi_EWpval,[],2);
b=reshape(dataOut.HGppEW,[],2);
subplot(2,2,3)
plot(dataOut.Np,-log10(a(TNidx,1)),'ko');hold on
plot(dataOut.Np,-log10(b(TNidx,1)),'r+');hold on
xlabel('Number of ROI pairs within network pair')
ylabel('FPR')
title(params.group{1})
[r,p]=corr(dataOut.Np,a(TNidx,1));
text(700,4,['CS r = ',num2str(r)])
text(700,3.5,['CS p = ',num2str(p)])
[r,p]=corr(dataOut.Np,b(TNidx,1));
text(700,3,['HG r = ',num2str(r)])
text(700,2.5,['HG p = ',num2str(p)])
ylim([0,5])
% axis([0,1000,1e-4,1])

subplot(2,2,4)
plot(dataOut.Np,-log10(a(TNidx,2)),'ko');hold on
plot(dataOut.Np,-log10(b(TNidx,2)),'r+');hold on
xlabel('Number of ROI pairs within network pair')
ylabel('FPR')
title(params.group{2})
[r,p]=corr(dataOut.Np,a(TNidx,2));
text(700,4,['CS r = ',num2str(r)])
text(700,3.5,['CS p = ',num2str(p)])
[r,p]=corr(dataOut.Np,b(TNidx,2));
text(700,3,['HG r = ',num2str(r)])
text(700,2.5,['HG p = ',num2str(p)])
ylim([0,5])
% axis([0,1000,1e-4,1])







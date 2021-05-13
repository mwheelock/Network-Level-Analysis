function Perm_2group_Ttest_1tp_Figs2(dataOut,IM,params)

% This function is the same as Perm_2group_Ttest_1tp_Figs2 except 
% that it does not plot the block-wise FPR values. 
% Also, it plots the FPR-p-vals for both
% Chi-squared and HG to see dynamic range

Pth=0.05/params.B;
PthMcN=0.05;
[NNidx,Tidx,TNidx]=IM2idx_HSB(IM);
Npairs=length(TNidx);

if (max(abs(dataOut.rho(:)))<=1)
    enMax=0.3;
else
    enMax=3;
end


Nnets=max(IM.key(:,2));
NetKey=[[1:Nnets]',[1:Nnets]'];
cMap=[1,1,1;0,0,0];
cMap2=[1,1,1;1,0,0;0,0,1;1,0,1];


%% Enrichment 
a=reshape(dataOut.Chi_pval0,[],1);b=reshape(dataOut.Chi_EWpval,[],1);
c=reshape(dataOut.HGp,[],1);d=reshape(dataOut.HGppEW,[],1);
for j=1
figure('Color','w','Position',[100,100,1400,800]);
subplot(2,3,1)
Matrix_Org3_HSB(squeeze(dataOut.rho(:,:,j)),IM.key,10,[-1,1].*enMax,IM.cMap,0,jet);
title([params.group{1} 'vs' params.group{2},'mean difference'])
subplot(2,3,4)
Matrix_Org3_HSB(squeeze(dataOut.Npval_lt_0p05(:,:,j)),IM.key,10,...
    [0,1],IM.cMap,0,cMap);
title([' p<',num2str(params.Pmax)])

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



switch params.scale
    case 'pvalue'
newmap = colormap(flipud(parula(1000)));
ncol = size(newmap,1);
newmap(ncol,:) = [1 1 1];
tempmap=colormap(newmap); 

subplot(2,3,5)
foo=(dataOut.Chi_EWpval);
Matrix_Org3_HSB(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx))-.001,Pth],IM.cMap,0,tempmap);
[r,col]=ind2sub(size(squeeze(dataOut.Chi_EWpval(:,:,j))),...
    find(squeeze(dataOut.Chi_EWpval(:,:,j))<Pth));
titlename = strcat('Chi-Square permuted p<',num2str(Pth));
title(titlename)
hold on
plot(col,r,'*k');colorbar

subplot(2,3,6)
foo=(dataOut.HGppEW);
Matrix_Org3_HSB(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx))-0.001,Pth],IM.cMap,0,tempmap);
[r,col]=ind2sub(size(squeeze(dataOut.HGppEW(:,:,j))),...
find(squeeze(dataOut.HGppEW(:,:,j))<Pth));
titlename = strcat('Hypergeometric permuted p<',num2str(Pth));
title(titlename)
hold on
plot(col,r,'*k');colorbar
   
        
    case 'neglog10'
subplot(2,3,5);
foo=-log10(dataOut.Chi_EWpval(:,:,j));
Matrix_Org3_HSB(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx)),2],IM.cMap,0,parula(1000));
[r,col]=ind2sub(size(squeeze(dataOut.Chi_EWpval(:,:,j))),...
    find(squeeze(dataOut.Chi_EWpval(:,:,j))<Pth));
hold on
title('Chi-Square scaled -log10(2) = p<0.01')
plot(col,r,'*k'); colorbar

        
subplot(2,3,6)
foo=-log10(dataOut.HGppEW(:,:,j));
Matrix_Org3_HSB(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx)),2],IM.cMap,0,parula(1000));
[r,col]=ind2sub(size(squeeze(dataOut.HGppEW(:,:,j))),...
    find(squeeze(dataOut.HGppEW(:,:,j))<Pth));
hold on
title('Hypergeometric scaled -log10(2) = p<0.01')
plot(col,r,'*k'); colorbar
end


%% Significant stats by network block
cMap2=[1,1,1;0.5,0.75,1;0,0,1;...% none, EN12(noMnC),EN12(McN)
    1,0.5,0;1,0,0;... %EN24(noMcN),EN24(McN)
    0,1,0;1,0,1]; %EN12,24(noMcN),EN12,24(McN)
Xpos=8;
FS=16;
Tcol='k';
Fcol='w';
E12=(squeeze(dataOut.Chi_EWpval(:,:,1))<=Pth).*...
    (squeeze(dataOut.HGppEW(:,:,1))<=Pth);
E24=zeros(size(dataOut.HGppEW));
McN=0<=Pth;
a=(E12.*(-E24+1).*(-McN+1)).*1+...% E12, no E24, no McN
    (E12.*(-E24+1).*(McN)).*2+...% E12, no E24,  McN;
    (E24.*(-E12+1).*(-McN+1)).*3+...% E24, no E12, no McN;
    (E24.*(-E12+1).*(McN)).*4+...% E24, no E12,  McN;
    (E12.*(E24).*(-McN+1)).*5+...% E12, E24, no McN;
    (E12.*(E24).*(McN)).*6; % E12, E24, McN;
figure('Color',Fcol,'Position',[100,100,900,600]);
Matrix_Org3_HSB(a,NetKey,0.5,[0,6],IM.cMap,0,cMap2);
hold on;
rectangle('Position',[Xpos,0.5,0.5,0.5],'FaceColor',cMap2(2,:))
text(Xpos+1,0.75,'Enrichment at 12 mo','Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,1.5,0.5,0.5],'FaceColor',cMap2(3,:))
text(Xpos+1,1.75,'Enrichment at 12 mo & McNemar','Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,2.5,0.5,0.5],'FaceColor',cMap2(4,:))
text(Xpos+1,2.75,'Enrichment at 24 mo','Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,3.5,0.5,0.5],'FaceColor',cMap2(5,:))
text(Xpos+1,3.75,'Enrichment at 24 mo & McNemar','Color',Tcol,'FontSize',FS)
rectangle('Position',[Xpos,4.5,0.5,0.5],'FaceColor',cMap2(6,:))
text(Xpos+1,4.75,'Enrichment at 12 & 24 mo','Color',Tcol','FontSize',FS)
rectangle('Position',[Xpos,5.5,0.5,0.5],'FaceColor',cMap2(7,:))
text(Xpos+1,5.75,'Enrichment at 12 & 24 mo & McNemar','Color','k','FontSize',FS)
% set(gcf,'Color','k')

%% Double check Network-pair size is not a confound
% a=reshape(dataOut.Chi_EWpval,[],2);
% b=reshape(dataOut.HGppEW,[],2);
a=reshape(dataOut.Chi_pval0,[],1);
b=reshape(dataOut.HGp,[],1);
figure;
subplot(1,2,1)
plot(dataOut.Np,-log10(a(TNidx,1)),'ko');hold on
plot(dataOut.Np,-log10(b(TNidx,1)),'r+');hold on
xlabel('Number of ROI pairs within network pair')
ylabel('Asymptotic p-value')
title('12 mo')
[r,p]=corr(dataOut.Np,a(TNidx,1));
text(700,7,['CS r = ',num2str(r)])
text(700,6,['CS p = ',num2str(p)])
[r,p]=corr(dataOut.Np,b(TNidx,1));
text(700,5,['HG r = ',num2str(r)])
text(700,4,['HG p = ',num2str(p)])
% axis([0,1000,1e-4,1])

% subplot(2,2,2)
% plot(dataOut.Np,-log10(a(TNidx,2)),'ko');hold on
% plot(dataOut.Np,-log10(b(TNidx,2)),'r+');hold on
% xlabel('Number of ROI pairs within network pair')
% ylabel('Asymptotic p-value')
% title('24 mo')
% [r,p]=corr(dataOut.Np,a(TNidx,2));
% text(700,5,['CS r = ',num2str(r)])
% text(700,4,['CS p = ',num2str(p)])
% [r,p]=corr(dataOut.Np,b(TNidx,2));
% text(700,3,['HG r = ',num2str(r)])
% text(700,2,['HG p = ',num2str(p)])
% % axis([0,1000,1e-4,1])

a=reshape(dataOut.Chi_EWpval,[],1);
b=reshape(dataOut.HGppEW,[],1);
subplot(1,2,2)
plot(dataOut.Np,-log10(a(TNidx,1)),'ko');hold on
plot(dataOut.Np,-log10(b(TNidx,1)),'r+');hold on
xlabel('Number of ROI pairs within network pair')
ylabel('FPR')
title('12 mo')
[r,p]=corr(dataOut.Np,a(TNidx,1));
text(700,4,['CS r = ',num2str(r)])
text(700,3.5,['CS p = ',num2str(p)])
[r,p]=corr(dataOut.Np,b(TNidx,1));
text(700,3,['HG r = ',num2str(r)])
text(700,2.5,['HG p = ',num2str(p)])
ylim([0,5])
% axis([0,1000,1e-4,1])

% subplot(2,2,4)
% plot(dataOut.Np,-log10(a(TNidx,2)),'ko');hold on
% plot(dataOut.Np,-log10(b(TNidx,2)),'r+');hold on
% xlabel('Number of ROI pairs within network pair')
% ylabel('FPR')
% title('24 mo')
% [r,p]=corr(dataOut.Np,a(TNidx,2));
% text(700,4,['CS r = ',num2str(r)])
% text(700,3.5,['CS p = ',num2str(p)])
% [r,p]=corr(dataOut.Np,b(TNidx,2));
% text(700,3,['HG r = ',num2str(r)])
% text(700,2.5,['HG p = ',num2str(p)])
% ylim([0,5])
% axis([0,1000,1e-4,1])

end





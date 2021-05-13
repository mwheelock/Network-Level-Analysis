function Perm_1tp_Figs2(dataOut,IM,B,params)


%% Set parameters
if ~exist('B','var'),B=1;end
Pth=0.05./B;
[NNidx,Tidx,TNidx]=IM2idx(IM);

Nnets=max(IM.key(:,2));
NetKey=[[1:Nnets]',[1:Nnets]'];
cMap=[1,1,1;0,0,0];

%% Enrichment 
a=reshape(dataOut.Chi_pval0,[],1);b=reshape(dataOut.Chi_EWpval,[],1);
c=reshape(dataOut.HGp,[],1);d=reshape(dataOut.HGppEW,[],1);

figure('Color','w','Position',[100,100,1400,800]);
subplot(2,3,1)
Matrix_Org3(squeeze(dataOut.rho),IM.key,10,[-0.3,0.3],IM.cMap,0);
switch params.type
    case 'Pearson'
title(['\rho'])
    case 'Spearman'
title(['\rho'])
end
subplot(2,3,4)
Matrix_Org3(squeeze(dataOut.thresholdedROIpairs),IM.key,10,[0,1],IM.cMap,0,cMap);
title(['p<0.05'])

subplot(2,3,2)
loglog(squeeze(dataOut.Emp_FDR_CS(:,1)),squeeze(dataOut.Emp_FDR_CS(:,2)),'k');
hold on;
loglog([1e-20,1],[Pth,Pth],'b','LineWidth',2)
loglog([squeeze(dataOut.Chi_EWth),squeeze(dataOut.Chi_EWth)],...
    [1e-5,1],'r','LineWidth',2)
loglog(squeeze(a(TNidx)),squeeze(b(TNidx)),'*k');
axis([1e-15,1,1e-3,1])
title('\chi^2 p-values');xlabel('Asymptotic');
ylabel('Permutation-based FPR')

subplot(2,3,5)
foo=-log10(dataOut.Chi_EWpval);
Matrix_Org3(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx)),2],IM.cMap,0,parula(1000));
[r,col]=ind2sub(size(dataOut.Chi_EWpval),find(dataOut.Chi_EWpval<Pth));
hold on
plot(col,r,'*k');%colorbar

subplot(2,3,3)
loglog(squeeze(dataOut.Emp_FDR_HG(:,1)),...
    squeeze(dataOut.Emp_FDR_HG(:,2)),'k');hold on;
loglog([1e-20,1],[Pth,Pth],'b','LineWidth',2)
loglog([squeeze(dataOut.HG_EWth),squeeze(dataOut.HG_EWth)],...
    [1e-5,1],'r','LineWidth',2)
loglog(squeeze(c(TNidx)),squeeze(d(TNidx)),'*k');
axis([1e-15,1,1e-3,1])
title('Hypergeometric p-values');xlabel('Asymptotic');
ylabel('Permutation-based FPR')

subplot(2,3,6)
foo=-log10(dataOut.HGppEW);
Matrix_Org3(squeeze(foo),...
    NetKey,0.5,[min(foo(TNidx)),2],IM.cMap,0,parula(1000));
[r,col]=ind2sub(size(dataOut.HGppEW),find(dataOut.HGppEW<Pth));
hold on
plot(col,r,'*k');%colorbar

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
title(params.group{1})
[r,p]=corr(dataOut.Np,a(TNidx,1));
text(700,7,['CS r = ',num2str(r)])
text(700,6,['CS p = ',num2str(p)])
[r,p]=corr(dataOut.Np,b(TNidx,1));
text(700,5,['HG r = ',num2str(r)])
text(700,4,['HG p = ',num2str(p)])
% axis([0,1000,1e-4,1])

a=reshape(dataOut.Chi_EWpval,[],1);
b=reshape(dataOut.HGppEW,[],1);
subplot(1,2,2)
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
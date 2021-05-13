function fcBx_PercFC_Scat2(pairsOut,IM)

% This function creates a scatter plot of the proportion of positive fcMRI 
% values against fcBx correlation from the output of 
% pairsOut=View_Clear_Brain_ROI_Corr2a.
% tp is the time point of interest (set to 1 if only 1)
% IM is the infomap structure.

%% Parameters and Initializations

cMapN=hsv(2800);cMapN=cMapN(1401:2400,:);
cMapP=cat(1,summer(500),flip(autumn(500),1));

%% Make figure;
set0=pairsOut;
rho=set0(:,3);
Pfc=set0(:,4);
Ndots=length(rho);
for j=1:Ndots
    switch sign(rho(j))
        case 1
            val=round(1000*Pfc(j));
            if ~val, val=1;end
            dataColor=cMapP(val,:);
        case -1
            val=round(1000*Pfc(j));
            if ~val, val=1;end
            dataColor=cMapN(val,:);            
    end
plot(Pfc(j),rho(j),'o','MarkerFaceColor',dataColor,...
    'MarkerEdgeColor','k','MarkerSize',10);
hold on
end
ylabel(['Brain-behavior Correlation']);
if (max(abs(rho(:)))<=1)
xlabel('p(fc[z])>0');axis([0,1,-1,1])
else
xlabel('p(fc[z])>0');axis([0,1,-6,6])
end
[h,p,ci,s]=ttest(rho);
axis square

title([[{[IM.Nets{IM.key(set0(1,1),2)},'-',...
    IM.Nets{IM.key(set0(1,2),2)}]}],...])%,...
   [{['t(\rho)=',num2str(s.tstat),', p=',num2str(p)]}]])
plot([0,1],[0,0],'--k');
% plot([1/3,1/3],[-1,1],'--k');plot([2/3,2/3],[-1,1],'--k');
% Ntot=size(set0,1);
% text(1/6,0.8,num2str(sum((Pfc<1/3).*(rho>0))/(Ntot/6)))
% text(3/6,0.8,num2str(sum((Pfc>=1/3).*(Pfc<2/3).*(rho>0))/(Ntot/6)))
% text(5/6,0.8,num2str(sum((Pfc>=2/3).*(rho>0))/(Ntot/6)))
% text(1/6,-0.8,num2str(sum((Pfc<1/3).*(rho<0))/(Ntot/6)))
% text(3/6,-0.8,num2str(sum((Pfc>=1/3).*(Pfc<2/3).*(rho<0))/(Ntot/6)))
% text(5/6,-0.8,num2str(sum((Pfc>=2/3).*(rho<0))/(Ntot/6)))
function stats=ModularityHSB(clrs,rmat0,rth,params)
%
% This function calculates the modularity for a set of clustering
% assignments for a given rmatrix. This metric is based on Newman, 2004.
% It is assumed that the clrs per kden/rth are along each column on clrs.
% Ignore modules set to zero (unclassified).

%% Set parameters

[Nroi,Nkden]=size(clrs);
stats.modularity=zeros(Nkden,1);
rmat=zeros(size(rmat0),'single');
UDidx=triu(ones(Nroi),1)==1;
domod=1;

%% Calculate modularity
for j=1:Nkden
    rmat=+(rmat0>=rth(j)); % Binarize
    
    if domod==1
        %disp('Calculate Modularity');
        clusters=squeeze(clrs(:,j));
        Ugroups=unique(clusters);
%         if Ugroups(1)==0    % Unclassified == 0
%             Ugroups(1)=[];
%         else                % Unclassified == max
            switch params.usp
                case 'first'
                     Ugroups(1)=[];
                case 'last'
                     Ugroups(end)=[];
                case 'none'
                   Ugroups=Ugroups; % do nothing
            end
 %       end
    
        
        % Initialize
        g=size(Ugroups,1);
        e0=zeros(g,'single');
        
        % Calc strength of edges within and btwn groups
        for x=1:g
            a=(clusters==Ugroups(x));
            for y=1:g
                e0(x,y)=sum(sum(rmat(a,(clusters==Ugroups(y)))));
            end
        end
        
        % Calc modularity
        e0=e0/sum(e0(:));
        tr=trace(e0);
        e2=sum((sum(e0)).^2);
        stats.modularity(j)=tr-e2;
        
        %disp('Calculate Nedges');
        stats.Nedges(j)=sum(rmat(UDidx)~=0);
        %disp('Calculate kave');
        stats.kave(j)=(stats.Nedges(j)*2)/Nroi;
        %disp('Calculate Assortativity');
        degrees=sum(rmat)';
        [x,y]=find(triu(rmat,1));
        b=[degrees(x),degrees(y)]-1;
        c=corrcoef(b);
        stats.A(j)=c(1,2);
        
    else
        
        disp('Calculate Num Components etc (takes time)');
        M=diag(-sum(rmat))+rmat;    % go go gadget Ncomponents
        [~,D]=eig(M);
        D=abs(diag(D));
        stats.Nc(j)=sum(D(:)<1e-4); % Ncomponents
        D=zeros(Nroi);
        n=1;nPATH=rmat;L=(nPATH~=0);
        while find(L,1), D=D+n.*L;n=n+1;nPATH=nPATH*rmat;L=(nPATH~=0).*(D==0);end
        R=single(D~=0);
        stats.NnBc(j)=max(sum(R))/Nroi; % Nn in biggest component
        stats.Cdns(j)=sum(R(:))/(Nroi*Nroi); % Connectedness
    end
end

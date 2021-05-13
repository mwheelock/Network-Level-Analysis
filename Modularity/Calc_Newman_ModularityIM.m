function [Qn_bin, Qtb, QmscUsP, QmscnoUsP] = Calc_Newman_ModularityIM(fc,IM,params)
% This function calculates Newman's Q (2004) and Traag and Bruggeman Qtb (2009) based on a predefined network
% structure (Infomap solution saved in IM.key). 
% This code also calculates Newman's Q on unbinarized matrix using a
% function from MSCcodebase-Master (Gordon et al. 2017 Neuron)
%
% The default settings assume an edge density threshold of 2.5% for Newman's Q and fc data
% are z values. These settings can be modified using the params structure
% params options
% params.type = 'kden' or 'r'
% params.datatype = 'zmat' or 'rmat' 
% params.thr = [specify numerical value for threshold]
% params.usp = 'first' or 'last' or 'none' Are unassigned/unspecified ROI listed first or last in IM.nets
% 
% Qn_bin = Newman's Q calculated on binarized matrix with USp network removed
% QtB = Traag and Bruggeman Qtb (2009)
% QmscUsP = Newman's Q calculated using non-binarized matrix without removing UsP network (most similar to Evan's Neuron paper)
% QmscnoUsP = Newman's Q calculated using non-binarized matrix removing the UsP network
% UsP = Unspecified network
%
% MW 9-16-20
%%
if ~exist('params','var'), params=struct; end
if ~isfield(params,'type'), params.type='kden';end % alt 'r'
if ~isfield(params,'datatype'), params.datatype='zmat';end % alt 'rmat'
if ~isfield(params,'thr'), params.thr=0.025;end
if ~isfield(params,'usp'), params.usp='last';end 

switch params.datatype
    case 'zmat' % if fc data are z values
        rmat = FisherZ2R_HSB(fc(IM.order,IM.order,:));
        rmatAve = squeeze(mean(rmat,3));
    case 'rmat'
        rmat = fc(IM.order,IM.order,:);
        rmatAve=squeeze(mean(rmat,3));
end

% visualize matrix organization
 Matrix_Org3_HSB(rmatAve,IM.key,10,[-0.3,0.3],IM.cMap,1);

% Initialize outputs
[~,Nroi,~]=size(rmat);
NPE=Nroi*(Nroi-1)/2;
Nkden=length(params.thr);
UDidx=find(triu(ones(Nroi),1)==1);      % indices of unique conns
modularity=zeros(Nkden,1,'single');    % Modularity
rth=modularity;                         % r threshold
kdenth=modularity;                      % kden threshold

for j=1:Nkden

    rmat0=rmatAve;    
    
    % threshold matrix
    switch params.type
        case 'r'
            rmat0(rmat0<params.thr(j))=0;
            kdenth(j)=sum(rmat0(UDidx)~=0)/NPE;
            rth(j)=min(rmat0(UDidx));
        case 'kden'
            EL=ceil(params.thr(j)*NPE);
            rmat0=triu(rmat0,1);
            [v,idx]=sort(rmat0(:),'descend');
            v((EL+1):length(UDidx))=0;
            rmat0(idx)=v;
            rmat0=max(rmat0,rmat0'); %thresholded matrix
            rth(j)=v(EL);
            kdenth(j)=EL/NPE;
    end
    
    %remove Unnasigned ROI network
    nets = unique(IM.key(:,2));
           switch params.usp
                case 'first'
                     rmnet=nets(1); % remove net
                     x = IM.key(:,2) ~= rmnet;
                     rmatAvenoUSP=rmatAve(x,x);
                     rmat0noUSP=rmat0(x,x);
                case 'last'
                    rmnet=nets(end); % remove net
                    x = IM.key(:,2) ~= rmnet;
                    rmatAvenoUSP=rmatAve(x,x);
                    rmat0noUSP=rmat0(x,x);
               case 'none'
                    rmat0noUSP=rmat0;
           end
                
end
                    
stats=ModularityHSB(IM.key(:,2),rmatAve,rth,params); % ATE function for Newman 2004 removes Usp Network
Qn_bin = stats.modularity; 
% need to input thresholded matrix into MSCcodebase function %
[QnoUSP Qds Qds_Li C] = calc_modularity_TL(IM.key(x,2),rmat0noUSP); % remove Usp Network 
[QUSP Qds Qds_Li C] = calc_modularity_TL(IM.key(:,2),rmat0); % keep Usp Network
QmscnoUsP = QnoUSP;
QmscUsP = QUSP;
Qtb=ModularityIM(rmatAve,IM); % Traag and Bruggeman, 2009

end
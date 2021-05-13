function dataOut=fcBx_Enrich_1tp(fc,Bx,IM,params)

% Perform the correlation analysis (type = {Pearson, Spearman}), 
% and enrichment analyses using Chi-square and Hypergeometric stats.
% fc is the Fisher-z fc-data; Nseeds-Nseeds-Nsubj
% Bx is the behavior data; Nsubj
% IM is the ROI-to-module key.
% params.type determines the type of correlation {'Spearman','Pearson'}
% params.np sets number of iterations shuffling subject fc-Bx pairs
% params.Pmax sets p-cutoff for roi-roi to network block analyses
% 
% This is currently only optimized for Matlab 2015a or later


%% Parameters & Initialization
if ~isfield(params,'type'), params.type='Spearman';end
if ~isfield(params,'nnPmax'), params.nnPmax=0.05;end
if ~isfield(params,'Pmax'), params.Pmax=0.05;end
if ~isfield(params,'B'), params.B=1;end
type=params.type;
Ntests=params.np;
Pmax=params.Pmax; % ROI-pair level p-value
nnPmax=params.nnPmax; % network-level p-value
B=params.B;

% General Bookkeeping
[NNidx,Tidx,TNidx]=IM2idx(IM);
Nroi_pairs=size(Tidx,1);
Nnets=max(IM.key(:,2));
Npairs=size(TNidx,1);
Np=cell2mat(cellfun(@length,NNidx,'UniformOutput',0));
[~,idx]=sort(Np);
cols=jet(Npairs);
cols(idx,:)=cols;

% roi data
[~,Nroi,Nsubjs]=size(fc);
rho=zeros(Nroi_pairs,1);
pval0i=zeros(Nroi_pairs,1);

% Network data
Chi_statsP=ones(Npairs,Ntests,'single');  % Chi stats for Perms
HGP=ones(Npairs,Ntests,'single');           % HG for Perms

% Initialize Outputs
dataOut.rho=zeros(Nroi*Nroi,1,'single');
dataOut.pval=zeros(Nroi*Nroi,1,'single');
dataOut.thresholdedROIpairs=zeros(Nroi*Nroi,1,'single');
dataOut.Chi_stats=zeros(Nnets*Nnets,1,'single');  % Output Chi stats
dataOut.Chi_pval0=ones(Nnets*Nnets,1,'single');   % Output naive pval
dataOut.Chi_pval=ones(Nnets*Nnets,1,'single');        % Output pval
dataOut.Chi_EWpval=ones(Nnets*Nnets,1,'single');        % Output pval
dataOut.Chi_pval_gte=zeros(Nnets*Nnets,1,'single');   % Output GTE
dataOut.HGp=ones(Nnets*Nnets,1,'single');% Output HG p-value
dataOut.HGpp=ones(Nnets*Nnets,1,'single');% Output HG perm p-value
dataOut.HGppEW=ones(Nnets*Nnets,1,'single');% Output HG perm p-value
dataOut.Np=Np;

c=logspace(-20,0,10001)';
       
fc=reshape(fc,[],Nsubjs);
fc=fc(Tidx,:);


%% Run initial Chi-squared test
switch type    
    case 'Spearman'
        [rho,pval0i]=MySpearman(Bx,fc');
    case 'Pearson'
        [rho,pval0i]=corr(Bx,fc','type','Pearson');
end

dataOut.rho(Tidx)=rho;
dataOut.pval(Tidx)=pval0i;

% Chi-Squared stat(Tidx)
pval0=+(pval0i<=Pmax);
N_exp_sig_ratio_tot=sum(pval0)./Nroi_pairs;
As=cell2mat(cellfun(@(x) sum(pval0(x)),NNidx,'UniformOutput',0));
Es=(N_exp_sig_ratio_tot*Np);
Chi_stats=((As-Es).^2).*((Es.^-1)+((Np-Es).^-1));
chi_pval_gte=+(As>Es);                 % GT test
Chi_stats(~isfinite(Chi_stats))=0;
dataOut.thresholdedROIpairs(Tidx)=pval0;
dataOut.Chi_stats(TNidx)=Chi_stats;  
dataOut.Chi_pval_gte(TNidx)=chi_pval_gte;   

% Hypergeometric stat
dataOut.HGp(TNidx)=hygecdf(As,Nroi_pairs,sum(As),Np,'upper');  
clear Tidx pval0 pval0i rho


%% Set p-values
if Ntests==1 % get naive p-values and mult-comp tests for 1-dof-chi
dataOut.Chi_pval(TNidx)=normcdf(sqrt(Chi_stats),'upper').^chi_pval_gte;


else % Run Permutation test
    disp(['Running ',num2str(Ntests),' Permutations'])
    
    Chi_pval0=erfc(sqrt(Chi_stats)./sqrt(2)).^chi_pval_gte;  
    dataOut.Chi_pval0(TNidx)=Chi_pval0;

    disp(['Running parfor'])
    parfor j=1:Ntests
            switch type
                case 'Spearman'
                    [~,pval0p]=MySpearman(Bx(randperm(Nsubjs)),fc');
                case 'Pearson'
                    [~,pval0p]=corr(Bx(randperm(Nsubjs)),fc');
            end

            pval=+(pval0p<=Pmax);
            N_exp_sig_ratio_tot=sum(pval)/Nroi_pairs;
            As=cell2mat(cellfun(@(x) sum(pval(x)),NNidx,'UniformOutput',0));
            Es=(N_exp_sig_ratio_tot*Np);
            Chi_statsP(:,j)=normcdf(...
                sqrt(((((As-Es).^2).*((Es.^-1)+((Np-Es).^-1))))),...
                'upper').^(As>Es);
            HGP(:,j)=(hygecdf(As,Nroi_pairs,sum(As),Np,'upper')).^(As>Es);
    end
    
    % Clear out some stuff
    clear fc Bx Np As Es pval NNidx 
    
    % Label True Chi2 stat Permed-p-value: tests-->(PDF-->)CDF-->p-dist
    disp('Calculating Permutation-based p-values')
    Chi_statsP(~isfinite(Chi_statsP))=1; % IFF As==Es, stat=0,p=1.
    
    disp('Plotting Network-Block-wise FPR')
    figure('Color','w','Position',[50,50,1220,575])
    for j=1:Npairs
        % Chi-squared statistics
        [~,idx]=sort([Chi_pval0(j),(squeeze(Chi_statsP(j,:)))]);
        dataOut.Chi_pval(TNidx(j))=find(idx==1)/(1+Ntests);
        % Experiment Wide
        temp=squeeze(Chi_statsP(:));
        [~,idx]=sort([Chi_pval0(j);(temp(:))]);
        dataOut.Chi_EWpval(TNidx(j))=find(idx==1)/(1+Ntests*Npairs);
        
        % Hypergeometric p-values
        [~,idx]=sort([dataOut.HGp(TNidx(j)),(squeeze(HGP(j,:)))]);
        dataOut.HGpp(TNidx(j))=find(idx==1)/(1+Ntests);
        % Experiment Wide
        temp=squeeze(HGP(:));
        [~,idx]=sort([dataOut.HGp(TNidx(j));(temp(:))]);
        dataOut.HGppEW(TNidx(j))=find(idx==1)/(1+Ntests*Npairs);
        
        subplot(1,2,1)
        h=hist((squeeze(Chi_statsP(j,:))),c)';
        loglog(c,cumsum(h./sum(h(:))),'-','Color',cols(j,:));hold on
        subplot(1,2,2)
        h=hist((squeeze(HGP(j,:))),c)';
        loglog(c,cumsum(h./sum(h(:))),'-','Color',cols(j,:));hold on
    end
    
    % Empirical FDR: how often does any permed NN have a stat>Xobs?
    Chi_statsP=reshape(Chi_statsP,[],1);
    h=hist(Chi_statsP,c)';
    dataOut.Emp_FDR_CS=[c,cumsum(h./sum(h(:)))];
    subplot(1,2,1)
    loglog(squeeze(dataOut.Emp_FDR_CS(:,1)),...
        squeeze(dataOut.Emp_FDR_CS(:,2)),'k','LineWidth',2);
    axis square;axis([1e-20,1e0,1e-6,1e0])
    xlabel('Chi-squared asymptotic p-val');
    ylabel('Permutation-based FPR')
    hold on
    
    HGP=reshape(HGP,[],1);
    h=hist(HGP,c)';
    dataOut.Emp_FDR_HG=[c,cumsum(h./sum(h(:)))];    
    subplot(1,2,2)
    loglog(squeeze(dataOut.Emp_FDR_HG(:,1)),...
        squeeze(dataOut.Emp_FDR_HG(:,2)),'k','LineWidth',2);
    axis square;axis([1e-20,1e0,1e-6,1e0])
    xlabel('Hypergeometric asymptotic p-val');
    ylabel('Permutation-based FPR')
    hold on
    
    subplot(1,2,1)
    for j=1:Npairs
        loglog(Chi_pval0(j),dataOut.Chi_pval(TNidx(j)),...
            '*','Color',cols(j,:));
    end
    subplot(1,2,2)
    for j=1:Npairs
        loglog(dataOut.HGp(TNidx(j)),dataOut.HGpp(TNidx(j)),...
            'o','Color',cols(j,:));
    end
end

% Set Exp-Wide threshold value, make sure it is below nnPmax!
[~,idx1]=min(abs(nnPmax-squeeze(dataOut.Emp_FDR_CS(:,2))));
dataOut.Chi_EWth=dataOut.Emp_FDR_CS(idx1,1);
if squeeze(dataOut.Emp_FDR_CS(idx1,2))>nnPmax
    dataOut.Chi_EWth=dataOut.Emp_FDR_CS((idx1-1),1);
end
[~,idx2]=min(abs(nnPmax-squeeze(dataOut.Emp_FDR_HG(:,2))));
dataOut.HG_EWth=dataOut.Emp_FDR_HG(idx2,1);
if squeeze(dataOut.Emp_FDR_HG(idx2,2))>nnPmax
    dataOut.HG_EWth=dataOut.Emp_FDR_HG((idx2-1),1);
end

% Set binary significance matrices
dataOut.Chi_sig=(dataOut.Chi_pval0<=dataOut.Chi_EWth);%.*...
%     (dataOut.Chi_pval<=nnPmax);
dataOut.HG_sig=(dataOut.HGp<=dataOut.HG_EWth);%.*(dataOut.HGpp<=nnPmax);


% "OR"
dataOut.Enrich_sig=(dataOut.Chi_sig+dataOut.HG_sig)>0;
% "AND"
dataOut.Enrich_sig2=(dataOut.Chi_sig.*dataOut.HG_sig)>0;


%% Fix output shapes
dataOut.rho=reshape(dataOut.rho,Nroi,Nroi,[]);
dataOut.pval=reshape(dataOut.pval,Nroi,Nroi,[]);
dataOut.ROIpairthreshold=Pmax;
dataOut.thresholdedROIpairs=reshape(dataOut.thresholdedROIpairs,Nroi,Nroi,[]);
dataOut.Chi_stats=reshape(dataOut.Chi_stats,Nnets,Nnets,[]);
dataOut.Chi_pval=reshape(dataOut.Chi_pval,Nnets,Nnets,[]);
dataOut.Chi_EWpval=reshape(dataOut.Chi_EWpval,Nnets,Nnets,[]);
dataOut.Chi_pval0=reshape(dataOut.Chi_pval0,Nnets,Nnets,[]);
dataOut.Chi_pval_gte=reshape(dataOut.Chi_pval_gte,Nnets,Nnets,[]);
dataOut.HGp=reshape(dataOut.HGp,Nnets,Nnets,[]);
dataOut.HGpp=reshape(dataOut.HGpp,Nnets,Nnets,[]);
dataOut.HGppEW=reshape(dataOut.HGppEW,Nnets,Nnets,[]);
dataOut.Chi_sig=reshape(dataOut.Chi_sig,Nnets,Nnets,[]);
dataOut.HG_sig=reshape(dataOut.HG_sig,Nnets,Nnets,[]);
dataOut.Enrich_sig=reshape(dataOut.Enrich_sig,Nnets,Nnets,[]);
dataOut.Enrich_sig2=reshape(dataOut.Enrich_sig2,Nnets,Nnets,[]);

%% Other figures
Perm_1tp_Figs2(dataOut,IM,B,params);
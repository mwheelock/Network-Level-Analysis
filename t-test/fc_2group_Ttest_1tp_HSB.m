function dataOut=fc_2group_Ttest_1tp_HSB(fc,IM,params)

% Performs the Welch-Satterthwaite-corrected t-test
% and enrichment analyses using Chi-square and Hypergeometric stats.
% fc is the Fisher-z fc-data; Nseeds-Nseeds-Nsubj
% IM is the ROI-to-module key.
% params.type determines the type of correlation {'Spearman','Pearson'}
% params.np sets number of iterations shuffling subject fc-Bx pairs
% params.Pmax sets p-cutoff for roi-roi to network block analyses
%
% This function runs the permutation analysis for multiple fc enrichment
% analyses in parallel so that they maintain their relevant correlation 
% structure. Each set of initial data is contained in a cell of fc.
%
% fc{1} contains fcMRI data for group 1
% fc{2} contains fcMRI data for group 2
%
% The primary analysis focusses on Enrichment Analysis for the full set of
% subjects at time point 1 and (separately) time point 2, as well as the
% McNemar statistic between those sets of data. Permutation-based
% significance is correlated separately for those 3 cases.
%
% While output data is structured into 2 levels (Enrichment for time points
% 1 and 2 and then a separate set of arrays for the McNemar stat between 
% them), the 4 data sets are maintained for randomization.
% 
% This is currently only optimized for Matlab 2015a or later


%% Set up parameters & initialize outputs
if ~isfield(params,'type'), params.type='Spearman';end
if ~isfield(params,'nnPmax'), params.nnPmax=0.05;end
if ~isfield(params,'Pmax'), params.Pmax=0.05;end
Ntests=params.np;
Pmax=params.Pmax;
nnPmax=params.nnPmax;
Nsets=1;
Ngroups=2;
Nbx=zeros(Nsets,1);

for j=1:Ngroups  % number of subjs for each group
    Nbx(j)=size(fc{j},3);
end

% General bookkeeping
[NNidx,Tidx,TNidx]=IM2idx_HSB(IM);
Nnets=max(IM.key(:,2));
Nroi_pairs=size(Tidx,1);
Npairs=size(TNidx,1);
Np=cell2mat(cellfun(@length,NNidx,'UniformOutput',0));
dataOut.Np=Np;
clear key

% roi data
[~,Nroi,~]=size(fc{1});
% rho=zeros(Nroi_pairs,Nsets);
% pval0i=zeros(Nroi_pairs,1);

% Network data
Chi_stats=cell(Nsets,1);
chi_pval_gte=cell(Nsets,1);Chi_pval0=cell(1,Nsets);  % Chi stats for Perms
Chi_statsP=ones(Npairs,Ntests,Nsets,'single');  % Chi stats for Perms
HGP=ones(Npairs,Ntests,Nsets,'single');           % HG for Perms
% McNs=zeros(Nroi_pairs,Nsets);
% McNp=ones(Npairs,Ntests,'single');               % McNemar Perm vals


% Initialize output data
dataOut.rho=zeros(Nroi*Nroi,Nsets,'single');
dataOut.pval=zeros(Nroi*Nroi,Nsets,'single');
dataOut.Npval_lt_0p05=zeros(Nroi*Nroi,Nsets,'single');
dataOut.Chi_stats=zeros(Nnets*Nnets,Nsets,'single');  % Output Chi stats
dataOut.Chi_pval0=ones(Nnets*Nnets,Nsets,'single');   % Output naive pval
dataOut.Chi_pval=ones(Nnets*Nnets,Nsets,'single');        % Output pval
dataOut.Chi_EWpval=ones(Nnets*Nnets,Nsets,'single');        % Output pval
dataOut.Chi_pval_gte=zeros(Nnets*Nnets,Nsets,'single');   % Output GTE
dataOut.HGp=ones(Nnets*Nnets,Nsets,'single');% Output HG p-value
dataOut.HGpp=ones(Nnets*Nnets,Nsets,'single');% Output HG perm p-value
dataOut.HGppEW=ones(Nnets*Nnets,Nsets,'single');% Output HG perm p-value
% dataOut.McNemar_b=ones(Nnets*Nnets,1,'single');
% dataOut.McNemar_c=ones(Nnets*Nnets,1,'single');
% dataOut.McNemar_P=ones(Nnets*Nnets,1,'single');
% dataOut.McNemar_PP=ones(Nnets*Nnets,1,'single');
% dataOut.McNemar_PPEW=ones(Nnets*Nnets,1,'single');
dataOut.TNidx=TNidx;

c=logspace(-20,0,10001)';    
    
if ~exist('Pmax','var'),Pmax=0.05;end
[~,idx]=sort(Np);
cols=jet(Npairs);
cols(idx,:)=cols;

% Prepare data: remove duplicate data and organize big Bx dataset
for j=1:Ngroups
fc{j}=reshape(fc{j},[],Nbx(j));
fc{j}=fc{j}(Tidx,:);
end
fcAll=cat(2,fc{1},fc{2});


%% Run initial Chi-squared test
for j=1:Nsets
[t,pval0i]=TtestWelch(fc{1},fc{2});
dataOut.rho(Tidx,j)=t;
dataOut.pval(Tidx,j)=pval0i;

% Chi-Squared stat(Tidx)
pval0=+(pval0i<=Pmax);
N_exp_sig_ratio_tot=sum(pval0,1)./Nroi_pairs;
As=cell2mat(cellfun(@(x) sum(pval0(x)),NNidx,'UniformOutput',0));
Es=(N_exp_sig_ratio_tot*Np);
Chi_stats{j}=((As-Es).^2).*((Es.^-1)+((Np-Es).^-1));
chi_pval_gte{j}=+(As>Es);                 % GT test
Chi_stats{j}(~isfinite(Chi_stats{j}))=0;
dataOut.Npval_lt_0p05(Tidx,j)=pval0;
dataOut.Chi_stats(TNidx,j)=Chi_stats{j};  
dataOut.Chi_pval_gte(TNidx,j)=chi_pval_gte{j};   

% Hypergeometric stat
dataOut.HGp(TNidx,j)=hygecdf(As,Nroi_pairs,sum(As),Np,'upper');%
% McNs(:,j)=pval0;
end    

% McNemar stat
% [dataOut.McNemar_P(TNidx),dataOut.McNemar_b(TNidx),...
%     dataOut.McNemar_c(TNidx)]=McNemar_fast(McNs,NNidx);
% figure;plot(dataOut.McNemar_b(TNidx)./...
%     (dataOut.McNemar_b(TNidx)+dataOut.McNemar_c(TNidx)),...
%     -log10(dataOut.McNemar_P(TNidx)),'k*');
% xlabel('b/(b+c)');ylabel('-log_1_0(p)');title('McNemar')
clear Tidx pval0 pval0i rho



%% Set p-values
if Ntests==1 % get naive p-values and mult-comp tests for 1-dof-chi
for j=1:Nsets
dataOut.Chi_pval(TNidx,j)=normcdf(...
    sqrt(Chi_stats{j}),'upper').^chi_pval_gte{j};   % Alex T.
end

else % Run Permutation test
    disp(['Running ',num2str(Ntests),' Permutations'])
    
for j=1:Nsets
    Chi_pval0{j}=erfc(...
        sqrt(Chi_stats{j})./sqrt(2)).^chi_pval_gte{j};   % naive p-val
    dataOut.Chi_pval0(TNidx,j)=Chi_pval0{j};
end

    disp(['Running parfor'])
    parfor j=1:Ntests
        
        randL=randperm(sum(Nbx)); % Longitudinal randomization
%         McNs=zeros(Nroi_pairs,2);
        for G=1:Nsets
            
            [~,pval0p]=TtestWelch(fcAll(:,randL(1:Nbx(1))),fcAll(:,randL((1+Nbx(1)):end)));

            pval=+(pval0p<=Pmax);
            N_exp_sig_ratio_tot=sum(pval)/Nroi_pairs;
            As=cell2mat(cellfun(@(x) sum(pval(x)),NNidx,'UniformOutput',0));
            Es=(N_exp_sig_ratio_tot*Np);
            Chi_statsP(:,j,G)=normcdf(...
                sqrt(((((As-Es).^2).*((Es.^-1)+((Np-Es).^-1))))),...
                'upper').^(As>Es);
            HGP(:,j,G)=(hygecdf(As,Nroi_pairs,sum(As),Np,'upper')).^(As>Es);
%             McNs(:,G)=pval;
        end
        % McNemar stat
%         McNp(:,j)=McNemar_fast(McNs,NNidx);
    end
    
    % Clear out some stuff
    clear fc p* Bx Np As Es pval NNidx
    
    % Label True Chi2 stat Permed-p-value: tests-->(PDF-->)CDF-->p-dist
    disp('Calculating Permutation-based p-values')
    Chi_statsP(~isfinite(Chi_statsP))=1; % IFF As==Es, stat=0,p=1.
    
    
    
    disp('Plotting Network-Block-wise FDR')
    figure('Color','w','Position',[50,50,1220,575])
    for j=1:Npairs
        for k=1:Nsets
        % Chi-squared statistics
        [~,idx]=sort([Chi_pval0{k}(j),(squeeze(Chi_statsP(j,:,k)))]);
        dataOut.Chi_pval(TNidx(j),k)=find(idx==1)/(1+Ntests);
        % Experiment Wide
        temp=squeeze(Chi_statsP(:,:,k));
        [~,idx]=sort([Chi_pval0{k}(j);(temp(:))]);
        dataOut.Chi_EWpval(TNidx(j),k)=find(idx==1)/(1+Ntests*Npairs);
        
        % Hypergeometric p-values
        [~,idx]=sort([dataOut.HGp(TNidx(j),k),(squeeze(HGP(j,:,k)))]);
        dataOut.HGpp(TNidx(j),k)=find(idx==1)/(1+Ntests);
        % Experiment Wide
        temp=squeeze(HGP(:,:,k));
        [~,idx]=sort([dataOut.HGp(TNidx(j),k);(temp(:))]);
        dataOut.HGppEW(TNidx(j),k)=find(idx==1)/(1+Ntests*Npairs);
        
        
        subplot(2,3,k)
        h=hist((squeeze(Chi_statsP(j,:,k))),c)';
        loglog(c,cumsum(h./sum(h(:))),'-','Color',cols(j,:));hold on
        if k==1, title('Time Point 1')
        elseif k==2, title('Time Point 2')
        end
        subplot(2,3,3+k)
        h=hist((squeeze(HGP(j,:,k))),c)';
        loglog(c,cumsum(h./sum(h(:))),'-','Color',cols(j,:));hold on
        end
        
        % McNemar p-values
%         [~,idx]=sort([dataOut.McNemar_P(TNidx(j)),McNp(j,:)]);
%         dataOut.McNemar_PP(TNidx(j))=find(idx==1)/(1+Ntests);
%         [~,idx]=sort([dataOut.McNemar_P(TNidx(j));McNp(:)]);
%         dataOut.McNemar_PPEW(TNidx(j))=find(idx==1)/(1+Ntests*Npairs);
%         subplot(2,3,k+1)
%         h=hist((squeeze(McNp(j,:))),c)';
%         loglog(c,cumsum(h./sum(h(:))),'-','Color',cols(j,:));hold on
%         title('McMemar')
    end

    
    % Empirical FDR: how often does any permed NN have a stat>Xobs?
    Chi_statsP=reshape(Chi_statsP,[],Nsets);
    HGP=reshape(HGP,[],Nsets);
    
    for k=1:Nsets
        h=hist(Chi_statsP(:,k),c)';
        dataOut.Emp_FDR_CS(:,:,k)=[c,cumsum(h./sum(h(:)))];
        
        h=hist(HGP(:,k),c)';
        dataOut.Emp_FDR_HG(:,:,k)=[c,cumsum(h./sum(h(:)))];
        
        subplot(2,3,k)
        loglog(squeeze(dataOut.Emp_FDR_CS(:,1,k)),...
            squeeze(dataOut.Emp_FDR_CS(:,2,k)),'k','LineWidth',2);
        axis square;axis([1e-20,1e0,1e-6,1e0])
        xlabel('Chi-squared asymptotic p-val');
        ylabel('Permutation-based FPR')
        subplot(2,3,3+k)
        loglog(squeeze(dataOut.Emp_FDR_CS(:,1,k)),...
            squeeze(dataOut.Emp_FDR_HG(:,2,k)),'k','LineWidth',2);
        axis square;axis([1e-20,1e0,1e-6,1e0])
        xlabel('Hypergeometric asymptotic p-val');
        ylabel('Permutation-based FPR')
        
        subplot(2,3,k)
        for j=1:Npairs
            loglog(Chi_pval0{k}(j),dataOut.Chi_pval(TNidx(j),k),...
                '*','Color',cols(j,:));
        end
        subplot(2,3,3+k)
        for j=1:Npairs
            loglog(dataOut.HGp(TNidx(j),k),dataOut.HGpp(TNidx(j),k),...
                'o','Color',cols(j,:));
        end
       
    % Set Exp-Wide threshold value
    [~,idx1]=min(abs(nnPmax-squeeze(dataOut.Emp_FDR_CS(:,2,k))));
    dataOut.Chi_EWth(k)=dataOut.Emp_FDR_CS(idx1,1,k);
    if squeeze(dataOut.Emp_FDR_CS(idx1,2,k))>nnPmax
        dataOut.Chi_EWth(k)=dataOut.Emp_FDR_CS((idx1-1),1,k);
    end
    [~,idx2]=min(abs(nnPmax-squeeze(dataOut.Emp_FDR_HG(:,2,k))));
    dataOut.HG_EWth(k)=dataOut.Emp_FDR_HG(idx2,1,k);
    if squeeze(dataOut.Emp_FDR_HG(idx2,2,k))>nnPmax
        dataOut.HG_EWth(k)=dataOut.Emp_FDR_HG((idx2-1),1,k);
    end
    
    % Set binary significance matrices
    dataOut.Chi_sig(:,k)=(dataOut.Chi_pval0(:,k)<=dataOut.Chi_EWth(k));%.*...
%                             (dataOut.Chi_pval(:,k)<=nnPmax);
    dataOut.HG_sig(:,k)=(dataOut.HGp(:,k)<=dataOut.HG_EWth(k));%.*...;
%                             (dataOut.HGpp(:,k)<=nnPmax);
    end
    
    
    % McNemar
%     h=hist(McNp(:),c)';
%     dataOut.Emp_FDR_McN=[c,cumsum(h./sum(h(:)))];
%     subplot(2,3,k+1)
%     loglog(squeeze(dataOut.Emp_FDR_McN(:,1)),...
%         squeeze(dataOut.Emp_FDR_McN(:,2)),'k','LineWidth',2);
%     axis square;axis([1e-20,1e0,1e-6,1e0])
%     xlabel('McNemar asymptotic p-val');
%     ylabel('Permutation-based FPR')
%     for j=1:Npairs
%         loglog(dataOut.McNemar_P(TNidx(j)),dataOut.McNemar_PP(TNidx(j)),...
%             '*','Color',cols(j,:));
%     end
%     
%     [~,idx3]=min(abs(nnPmax-squeeze(dataOut.Emp_FDR_McN(:,2))));
%     dataOut.McN_EWth=dataOut.Emp_FDR_McN(idx3,1);
%     if squeeze(dataOut.Emp_FDR_McN(idx3,2))>nnPmax
%         dataOut.McN_EWth=dataOut.Emp_FDR_McN((idx3-1),1);
%     end
%     dataOut.McN_sig=(dataOut.McNemar_P<=dataOut.McN_EWth);%.*...;
%                             (dataOut.McNemar_PP<=nnPmax);    
    % OR
    dataOut.Enrich_sig=(dataOut.Chi_sig+dataOut.HG_sig)>0;
    % AND
    dataOut.Enrich_sig2=(dataOut.Chi_sig.*dataOut.HG_sig)>0;
end

%% Fix output shapes
dataOut.rho=reshape(dataOut.rho,Nroi,Nroi,[]);
dataOut.pval=reshape(dataOut.pval,Nroi,Nroi,[]);
dataOut.Npval_lt_0p05=reshape(dataOut.Npval_lt_0p05,Nroi,Nroi,[]);
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
% dataOut.McNemar_b=reshape(dataOut.McNemar_b,Nnets,Nnets);
% dataOut.McNemar_c=reshape(dataOut.McNemar_c,Nnets,Nnets);
% dataOut.McNemar_P=reshape(dataOut.McNemar_P,Nnets,Nnets);
% dataOut.McNemar_PP=reshape(dataOut.McNemar_PP,Nnets,Nnets);
% dataOut.McNemar_PPEW=reshape(dataOut.McNemar_PPEW,Nnets,Nnets);
% dataOut.McN_sig=reshape(dataOut.McN_sig,Nnets,Nnets);

%% Other figures
% Perm_Mixed_LCS_Figs2(dataOut,IM,1);
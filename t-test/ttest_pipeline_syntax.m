%% Run t-test Enrichment
TCfc= fc(:,:,(group==1));
VPTfc = fc(:,:,(group==2));

fccombined={ TCfc VPTfc};
Bxname='group1 vs group2';
group1 = 'group1';
group2 = 'group2';

% Set parameters
params.np=1e4;        % Number of permutations: set to 1 to get simple stats
%params.type='Pearson'; % flavor of correlations: 'Pearson' or 'Spearman'
params.Pmax=0.05;
params.nnPmax=0.05; 
params.group = {group1 group2};
params.BxName=Bxname;
params.B=1; % Bonferonni level (i.e., how many Bx are you using?)
params.fn=strcat(Bxname,'_',num2str(params.np),'_',num2str(params.Pmax),'_',IM.name,'.mat'); 
%params.scale='pvalue'; 

% no covariates %
dataOut=fc_2group_Ttest_1tp_HSB(fccombined,IM,params)
Perm_2group_Ttest_1tp_Figs2(dataOut,IM,params)
save(params.fn,'dataOut','fc','params');


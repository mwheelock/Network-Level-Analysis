function  [rrmat, stats, QCFC, QCFCed]=plot_distance_QC(fc,motion,distances);
%% Plot correlation coefficient-QC vs edge distance 

Nsubjs=size(fc,3);
Nroi=size(fc,1);
rr=zeros(Nroi*Nroi,1);
connectivitydata=reshape(fc,[],Nsubjs);%fc-X-subjs
rr=corr(connectivitydata',motion,'type','Pearson');
[rr p]=corr(connectivitydata',motion,'type','Pearson');
rrmat = reshape(rr,[Nroi,Nroi]);
upper = triu(rrmat); upper(isnan(upper))=[];upper(upper==0)=[];
QCFC = median(abs(upper));
 %% Plot the data and make meaningful? statistics
nodiagrrmat = upper;
nodiagrrmat(isnan(nodiagrrmat))=[]; nodiagrrmat(nodiagrrmat==0)=[];
scatter(distances,nodiagrrmat);
lsline
[QCFCed pd] = corrcoef(distances,nodiagrrmat);
stats = regstats(distances,nodiagrrmat, 'linear');

end

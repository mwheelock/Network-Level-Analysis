function [fcr, Bxr] = partialvariance(fc,Bx,covariates)
%% Control for Covariates %%
% covariates are cov1, cov2... covN assumed to be column vectors
%MW 1-3-2020
%%
Nroi=size(fc,1);
Nsubs = size(fc,3);
sizecov = size(covariates,2);

meancovariates = zeros(Nsubs,sizecov);
% mean center covariates
for i=1:sizecov
cov = covariates(:,i); % raw
covmean = cov-mean(cov); % mean center
meancovariates(:,i) = covmean;
end

% mean center the fc 
fctemp = reshape(fc,[],Nsubs);
%meanfc = bsxfun(@minus,fctemp,mean(fctemp,2)); %mean centering fc makes it difficult to interpret - suggest not doing this
meanfc = fctemp;

% Estimate the residual fc partialling out covariates
Betafc = pinv(meancovariates) * meanfc';
fcrtemp = (meanfc'-meancovariates * Betafc)'; %fcr is the residual fc
fcr = reshape(fcrtemp,Nroi,Nroi,Nsubs);

% Partial out covariates from Bx of interest
meanBx = Bx-mean(Bx); %mean center
BetaBx = pinv(meancovariates)*meanBx;
Bxr = meanBx-meancovariates*BetaBx; % Bxr is the residual Bx

end
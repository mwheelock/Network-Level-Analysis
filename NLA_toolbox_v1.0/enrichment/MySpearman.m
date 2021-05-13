function [rho,pval]=MySpearman(x,y)

% This function is to be a dedicated super-fast Spearman correlation
% calculation for the input data X(col vector) and Y(matrix). Inputs are
% assumed to have the same number of rows and have no NaNs or Infs. Y is
% also assumed to have no duplicates within a column.
% Steps:    
%   1.  Sort data sets: X-->x, Y-->y. Assume ties exist.
%   2.  Calculate Pearson correlation coefficient of rankings = rho.
%       a. subtract mean
%       b. norm
%       c. dot product
%   3.  Calculation of p-value from either
%       a. Fisher value: z=sqrt((n-3)/(1.06))*arctanh(rho);
%       b. Student's t:  t=rho*sqrt((n-2)/(1-rho^2)); dof=n-2;

%% Parameters
n=size(x,1);%,p1
p2=size(y,2);
n3const = (n+1)*n*(n-1) ./ 3;
% rho = zeros(p1,p2,'single');
% pval = zeros(p1,p2,'single');    
        
%% Sort
[x,xadj] = tiedrank(x,0);
% [y1,yadj] = tiedrank(y,0);
[~,idx]=sort(y);
for j=1:p2, y(idx(:,j),j)=[1:n]';end

%% Calc rho and p-val
% for i = 1:p1    % for each column in X
D = sum(bsxfun(@minus,x,y).^2); % sum((xranki - yrankj).^2);

% meanD = (n3const - (xadj+yadj)./3) ./ 2;
% stdD = sqrt((n3const./2 - xadj./3)*(n3const./2 - yadj./3)./(n-1));
% % ASSUMING yadj==0!
meanD = (n3const - (xadj)./3) ./ 2;
stdD = sqrt((n3const./2 - xadj./3)*(n3const./2)./(n-1));

% rho
rho = ((meanD - D) ./ (sqrt(n-1)*stdD))'; 
% Limit off-diag correlations to [-1,1].
% t = find(abs(rho) > 1); 
% rho(t) = rho(t)./abs(rho(t)); % preserves NaNs

% p-val: Use a t approximation.
t = Inf*sign(rho);
ok = (abs(rho) < 1);
t(ok) = rho(ok) .* sqrt((n-2)./(1-rho(ok).^2));
pval = (2*tcdf(-abs(t),n-2));
% end
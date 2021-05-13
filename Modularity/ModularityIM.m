function Q=ModularityIM(A,IM)
%
% This function calculates the modularity of a matrix that is assumed to be
% Pearson correlation values. The structure of the matrix is described by
% the IM structure. The matrix, A, is assumed to be square, 
% have both positive and negative values and be bound by +/-1. 
% Calculations are following Traag and Bruggeman, 2009.
%   m == edges --> correlations
%   m_ors=m_ors_p+m_ors_n --> split into pos/neg parts (obs,exp)
%   Q{sigma}=(1/m)*sum[m_oss-m_ess]|communities (don't double count edges!)
%   A=A_p+A_n
%   m_ors_p=sum[A_p]|within community-pair r,s (m_ors_n w A_n)
%   m_ers_p=sum[Prs_p]|within community-pair r,s (m_ers_n w Prs_n)
%   Pij_p=k_p/m_p_t --> pos edges within networks over total num pos
%                           edges
% H{sigma}=sum{a_ss}|over all communities s
%         =-sum{a_rs}|over all communite pairs such that r~=s
%   a_rs=(m_ors_p-m_ors_n)-(m_ers_p-m_ers_n); a == Adhesion
% Roughly speaking, Q~-(1/m_tot)H.
%
% Alternatively, 
% H{sigma}=-sum[Aij-(g_p*Pij_p-g_n*Pij_n)]|over i==j (within community)

%% Parameters and Initialization
g_p=1;
g_n=1;

Nroi=size(A,1);
Nnets=max(IM.key(:,2));
Tidx=find(tril(ones(Nroi),-1)==1);
TNidx=find(tril(ones(Nnets),0)==1);
Npairs=size(TNidx,1);
rs=repmat(squeeze(IM.key(:,2)),[1,Nroi]);
cs=rs';
NNTidx=[rs(Tidx),cs(Tidx)];
clear rs cs
NNidx=cell(Npairs,1); % Each cell is a network pair assoc to Tidx 
n=0;
Win=[];
Wout=[];
for j=1:Nnets
    for k=j:Nnets
        n=n+1;
        if j==k, Win=cat(1,Win,n);else Wout=cat(1,Wout,n);end
        NNidx{n,1}=intersect(find(NNTidx(:,1)==k),find(NNTidx(:,2)==j));
    end
end
Np=cell2mat(cellfun(@length,NNidx,'UniformOutput',0));


%% Modularity

% Random weight probs
A_p=A.*(A>0);
A_p(A_p==0)=NaN;
ki_p=nansum(A_p,1);
kj_p=nansum(A_p,2);
A_n=-A.*(A<0);
A_n(A_n==0)=NaN;
ki_n=nansum(A_n,1);
kj_n=nansum(A_n,2);
m_p_t=nansum(A_p(:))/2;    % Total pos edges mag
m_n_t=nansum(A_n(:))/2;    % Total neg edges mag
Pij_p=bsxfun(@times,ki_p,kj_p)./(2*m_p_t);
Pij_n=bsxfun(@times,ki_n,kj_n)./(2*m_n_t);

% Prep for community-pair-wise calculations
Pij_p=Pij_p(Tidx);
Pij_n=Pij_n(Tidx);
A=A(Tidx);
% Pij_p=repmat(m_p_t/(Nroi*(Nroi-1)),size(A));
% Pij_n=repmat(m_n_t/(Nroi*(Nroi-1)),size(A));

% Hamiltonian 1
Harg=cell2mat(cellfun(@(x) sum((A(x)-(g_p.*Pij_p(x)-g_n*Pij_n(x)))),...
    NNidx,'UniformOutput',0));
H=-sum(Harg(Win));

% Modularity
Q=-(1/(m_p_t+m_n_t))*H;


% Expected
 m_ers_p=cell2mat(cellfun(@(x) nansum(A_p(x))./m_p_t,NNidx,'UniformOutput',0));
 m_ers_n=cell2mat(cellfun(@(x) nansum(A_n(x))./m_n_t,NNidx,'UniformOutput',0));
% 
% % Observed
 m_ors_p=cell2mat(cellfun(@(x) nansum(A_p(x)),NNidx,'UniformOutput',0));
 m_ors_n=cell2mat(cellfun(@(x) nansum(A_n(x)),NNidx,'UniformOutput',0));
% 
% % Adhesion
 a_rs=(m_ors_p-m_ors_n)-(m_ers_p-m_ers_n);

% H=sum(a_rs(Win));








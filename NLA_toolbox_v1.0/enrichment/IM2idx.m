function [NNidx,Tidx,TNidx]=IM2idx(IM)
%
% This function provides indexing into module pairs for fast reference based
% on the IM structure of a grouped network.

Nroi=length(IM.key);
Nnets=max(IM.key(:,2));

Tidx=find(tril(ones(Nroi),-1)==1);
TNidx=find(tril(ones(Nnets),0)==1);% Lower triangle w main diag for NN-pairs

Npairs=size(TNidx,1);
rs=repmat(IM.key(:,2),[1,Nroi]);
cs=rs';
NNTidx=[rs(Tidx),cs(Tidx)];
clear rs cs
NNidx=cell(Npairs,1); % Each cell is a network pair assoc to Tidx 
n=0;
for j=1:Nnets
    for k=j:Nnets
        n=n+1;
        NNidx{n,1}=intersect(find(NNTidx(:,1)==k),find(NNTidx(:,2)==j));
    end
end

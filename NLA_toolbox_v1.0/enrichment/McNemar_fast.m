function [p,b,c,CS]=McNemar_fast(d1,bk)
%
% This function is optimized to run McNemar on 2 sets of data that are
% already binarized. The bk array contains bookkeeping information to keep
% track of functional networks.

Npairs=length(bk);
for j=1:Npairs
    b(j)=sum((d1(bk{j,1},1)).*(~d1(bk{j,1},2)));
    c(j)=sum((~d1(bk{j,1},1)).*(d1(bk{j,1},2)));
end

bPc=b+c;
keep=bPc>25;

stat1=((b-c).^2)./bPc;
stat2=((abs(b-c)-1).^2)./bPc;

CS(keep)=stat1(keep);     % if b+c>25
CS(~keep)=stat2(~keep);   % if b+c<=25
CS(~isfinite(CS))=0;      % if both b and c are zero

p=normcdf(sqrt(CS),'upper');  
p(bPc==0)=1;
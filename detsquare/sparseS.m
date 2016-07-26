function [S]=sparseS(L,r)
S=zeros([L 1]);
coefs=rand([r 1]);
S(1:r,1)=coefs;
S(:,1)=S(randperm(L),1);

S=S./sum(S);
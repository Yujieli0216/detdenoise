function [Subij]=sub(S,i,j)

[n,m]=size(S);
A=zeros(n,m);

        Hij=S;
        Hij(i,:)=[];
        Hij(:,j)=[];
        A(i,j)=(-1)^(i+j)*det(Hij);
        Subij=Hij;
 
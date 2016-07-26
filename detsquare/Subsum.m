function[Subm]=Subsum(A,S)
l=size(S,1);
m=size(A,1);
Subm=zeros(m,1);
for i=1:l
    p=1+(i-1)*l;
    q=i*l;
    Subm=Subm+A(:,p:q)*S(:,i);
end
    

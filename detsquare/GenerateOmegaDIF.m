function [Omega] = GenerateOmegaDIF(n)
Omega=zeros(n^2,n^2);
count=1;
for k=1:1:n
    for j=1:1:n-1,
        Image=zeros(n,n);
        Image(k,j)=1; 
        Image(k,j+1)=3; 
        Omega(count,:)=Image(:)';
        count=count+1;
    end
%     Image=zeros(n,n);
%     Image(k,n)=1; 
%     Image(k,1)=0; 
%     Omega(count,:)=Image(:)';
%     count=count+1;
% end
% for k=1:1:n
%     for j=1:1:n-1,
%         Image=zeros(n,n);
%         Image(j,k)=2; 
%         Image(j+1,k)=0; 
%         Omega(count,:)=Image(:)';
%         count=count+1;
%     end
    Image=zeros(n,n);
    Image(n,k)=2; 
    Image(1,k)=1; 
    Omega(count,:)=Image(:)';
    count=count+1;
end
Omega=normcols(Omega);
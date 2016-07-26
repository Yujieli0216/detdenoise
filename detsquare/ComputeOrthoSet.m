function [Q, inds] = ComputeOrthoSet(Q0)
Q=Q0;
cnt=1;
inds=zeros(1,size(Q,1));
inds(cnt)=1;

for j=2:size(Q,1)
    Q(j,:)=Q(j,:)-(Q(j,:)*Q(1:cnt,:)')*Q(1:cnt,:);
    v=sqrt(Q(j,:)*Q(j,:)');
    if v>1e-6
        cnt=cnt+1;
        Q(cnt,:)=Q(j,:)/v;
        inds(cnt)=j;
    end
end
Q=Q(1:cnt,:);
inds=inds(1:cnt);
        
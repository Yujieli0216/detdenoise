function [A_est, W, S_est, obj]=IVM_QP(X0,r)
% [1] Coleman, T.F. and Y. Li, "A Reflective Newton Method for Minimizing a Quadratic Function Subject to Bounds on some of the Variables," 
% SIAM Journal on Optimization, Vol. 6, Number 4, pp. 1040-1058, 1996.



% beta=max(max(X0));  %numerical processing 

% while beta>=1
X0=diag(1./sum(X0,2))*X0;
% % % beta0=max(max(X0));
% % % X0=X0./beta0;
[mm0 nn0]=size(X0);
if mm0<3
    [A_est S_est obj]=nLCA_IVM2(X0);
    return;
else
delt=0;
X1=X0';
% % % % % % [K vol]=convhulln(X1);
X=(red_const_qhull(X1))';
% X=X1;
% m=length(X0(:,1));
% [A_e, indice, Rp ]= VCA(X0,'Endmembers',r,'SNR',25,'verbose','off'); 
% X=A_e;


[m0 n0]=size(X);
XX=X0*X0';
m=mm0;
% XX=eye(m,m);

% W=eye(mm0,mm0);
rand('state',10);
W=randn(r,mm0);
% w=diag(1./sum(w,2))*w;
v=det(W*W');
j00=1;
k00=1;
iter=1;
XX0=XX
beta=1;
% % beta=1/max(max(XX0));
% beta=10000;
% kkkk
while iter<=4*r
% for iter=1:1*r
    iter;
    XX=beta*XX0;
    WW=W*XX*W';
    A_P=adj_mat(WW);
    i=mod(iter,r);
    if i==0 i=r; end
    coef11=zeros(m,m);
    coef21=zeros(m,m);
    coef31=zeros(m,m);
    for j=1:r
        temp1(:,:)=A_P(i,j,:,:);
        A_P_temp1=adj_mat(temp1);

        if j<i
            coef10=zeros(1,m);
            for t=1:r-1
                 temp2(:,:)=A_P_temp1(t,i-1,:,:);
                 tt=t;
                if t>=i tt=t+1; end
                coef10=coef10+(-1)^(t+i-1)*det(temp2)*W(tt,:)*XX;                
            end
            coef11=coef11+(-1)^(i+j)*XX*W(j,:)'*coef10;
        end
        if j==i
            coef21=XX*det(temp1);
        end
        if j>i
            coef30=zeros(1,m);
            for t=1:r-1
                 temp3(:,:)=A_P_temp1(t,i,:,:);
                 tt=t;
                if t>=i tt=t+1; end
                coef30=coef30+(-1)^(t+i)*det(temp3)*W(tt,:)*XX;                
            end
            coef31=coef31+(-1)^(i+j)*XX*W(j,:)'*coef30;
        end
        obj_coef=coef11+coef21+coef31;
    end

    Aeq=ones(1,m0);
    beq=1;
    H=(obj_coef+obj_coef')/2;
    f=zeros(m,1);
    A=-X';
    b=zeros(n0,1);
    % LB=zeros(m0,1);
%     [y1,fval1,exitflag1]=linprog(f,A01,b,Aeq,beq,[],[],W(i,:)',optimset('Display','off'));
    [y,fval,exitflag]=quadprog1(-H,f,A,b,Aeq,beq,[],[],W(i,:)',optimset('Display','off'));
%     i;

%     y2
%     fval2
%     exitflag1
%     exitflag2

%     p=f'*y2; 
    obj(1,iter)=exitflag; 
    obj(2,iter)=fval; 
    if exitflag<0
        beta=beta*10;
        iter=iter-1;        
    end
    

    if exitflag>-2
        W(i,:)=y';
%         if abs(p)>abs(q) W(i,:)=y2';
%         else W(i,:)=y1';
%         end
    end
% % %     p=f'*W(i,:)'; 
% % % %     if i==mm0 obj01(j00)=abs((v-max(abs(p),abs(q))))/v; j00=j00+1; v=max(abs(p),abs(q)); end
% % %     obj01(j00)=v; j00=j00+1; v=p; 
    y;
   W;
   iter=iter+1;
end
W
if r==length(X0(:,1))
    A_est=inv(W);
else 
    A_est=zeros(r,r);
end
  S_est=W*X0;

beta

    
    

% for i=1:2*mm0
%     
%     A=-X';
%     b=zeros(n0,1);
%     Aeq=ones(1,m0);
%     beq=1;
%     % LB=zeros(m0,1);
%     [y1,fval1,exitflag1]=linprog(f,A,b,Aeq,beq,[],[],W(kk,:)',optimset('Display','off'));
%     [y2,fval2,exitflag2]=linprog(-f,A,b,Aeq,beq,[],[],W(kk,:)',optimset('Display','off'));
% %     i;
% %     y2;
% %     fval2;
% %     exitflag2;
%     p=f'*y2; 
%     obj0(1,k0)=exitflag1; 
%     q=f'*y1; obj00(2,k0)=exitflag2; k0=k0+1;
%     if abs(p)>abs(q) W(kk,:)=y2';
%     else W(kk,:)=y1';
%     end
%     if kk==mm0 obj01(j)=abs((v-max(abs(p),abs(q))))/v; j=j+1; v=max(abs(p),abs(q)); end
% end

% figure
% plot(obj01)
figure
subplot(2,1,1)
plot(obj(1,:))
subplot(2,1,2)
plot(obj(2,:))
% A_est


end

function [A, S, obj]=nLCA_IVM2(X0)
obj=1;
X0=diag(1./sum(X0,2))*X0;
[mm0 nn0]=size(X0);
temp_X=X0(1,:)-X0(2,:);
ind1=find(temp_X>0);
ind2=find(temp_X<0);
X_ind1=X0(:,ind1);
X_ind2=X0(:,ind2);
temp1=-X_ind1(2,:)./temp_X(ind1);
temp2=-X_ind2(2,:)./temp_X(ind2);
W11=max(temp1);
W21=min(temp2);
W12=1-W11;
W22=1-W21;
WW=[W11 W12; W21 W22];
A=inv(WW);
% for i=1:2
%     ind=find(A(i,:)<0);
%     A(i,ind)=0;
% end
S=WW*X0;

function Xf=red_const_qhull(X)

[L,M]=size(X);
index=find(sum(abs(X'))>=1e-6);
Xn=X(index,:); 
LL=length(index);
Y=sum(Xn');
Xh=Xn./(Y'*ones(1,M));
K = convhulln(Xh,{'QJ'});
% K = convhulln(Xh,{'QJ','Pp'});
[mm nn]=size(K);
index_s=sort(reshape(K,mm*nn,1));
% index_s=sort(vec(K)); % sorting

[leng,t]=size(index_s);
p=zeros(leng,1);
p(1,1)=1;
l=1;
while (l<leng)
    if index_s(l+1)~=index_s(l)
       p(l+1,1)=1;
    end
    l=l+1;
end
index_f=index_s(p>0); % list the indices of data corresponding to possible constraints
Xf=Xh(index_f,:);

function A_P=adj_mat(A)
[m n]=size(A);
if m>n||m<n
    fprintf('the matrix is not square\n');
    return;
end
for i=1:m
    ind1=[1:i-1 i+1:m];
    for j=1:m
        ind2=[1:j-1 j+1:m];
        A_P(i,j,:,:)=A(ind1,ind2);
%         B(:,:)=A_P(i,j,:,:)
    end
end

function [A_est W S_est  ratio recoveryH DET err]=IVM_QPR2(X0,r,Omega,H0)
% [1] Coleman, T.F. and Y. Li, "A Reflective Newton Method for Minimizing a Quadratic Function Subject to Bounds on some of the Variables," 
% SIAM Journal on Optimization, Vol. 6, Number 4, pp. 1040-1058, 1996.



% beta=max(max(X0));  %numerical processing 

% while beta>=1
[mm0 nn0]=size(X0);
X1=X0(:,1:nn0/2);
X1=diag(1./sum(X1,2))*X1;
X2=X0(:,(nn0/2+1):(nn0*2/2));
X2=diag(1./sum(X2,2))*X2;


% % % beta0=max(max(X0));
% % % X0=X0./beta0;

if mm0<3
    [A_est S_est obj]=nLCA_IVM2(X0);
    return;
else
delt=0;
X=X1;



[m0 n0]=size(X1);
XX=X1*X1';
m=mm0;

rand('state',10);
W1=randn(r,mm0);
W2=randn(r,mm0);
% w=diag(1./sum(w,2))*w;
v=det(W1*W1');
j00=1;
k00=1;
iter1=1;
iter2=1;
XX0=XX;
beta=1;
aa=0.0000005;
DET1=zeros(3*r,1);
DET2=zeros(3*r,1);
DET3=zeros(3*r,1);
DET=zeros(3*r,1);
while iter1<3*r
% for iter=1:1*r
tic
    iter1;
    XX=beta*XX0;
    WW=W1*XX*W1';
    A_P=adj_mat(WW);
    i=mod(iter1,r);
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
                if t>=i 
                    tt=t+1; 
                end
                coef10=coef10+(-1)^(t+i-1)*det(temp2)*W1(tt,:)*XX;                
            end
            coef11=coef11+(-1)^(i+j)*XX*W1(j,:)'*coef10;
        end
        if j==i
            coef21=XX*det(temp1);
        end
        if j>i
            coef30=zeros(1,m);
            for t=1:r-1
                 temp3(:,:)=A_P_temp1(t,i,:,:);
                 tt=t;
                if t>=i 
                    tt=t+1;
                end
                coef30=coef30+(-1)^(t+i)*det(temp3)*W1(tt,:)*XX;                
            end
            coef31=coef31+(-1)^(i+j)*XX*W1(j,:)'*coef30;
        end
        obj_coef=coef11+coef21+coef31;
    end

    Aeq=ones(1,m0);
    beq=1;
    H=(obj_coef+obj_coef')/2;
    f=zeros(m,1);
    A=-X1';
    b=zeros(n0,1);
    
    [y,fval,exitflag]=quadprog1(-H,f,A,b,Aeq,beq,[],[],W1(i,:)',optimset('Display','off'));

    obj1(1,iter1)=exitflag; 
    obj1(2,iter1)=fval; 
    if exitflag<0
        beta=beta*1.2;
         iter1=iter1-1;        
    end
    

    if exitflag>-2
        W1(i,:)=y';
    end
    y;
    W1;
     iter1=iter1+1;
     A_est1=inv(W1);
      S_est1=W1*X1;
      if iter1~=0
        DET1(iter1)=-log(abs(det(S_est1*S_est1')));
 
      end    
      toc
end


XX=X2*X2';
m=mm0;

rand('state',10);

W2=randn(r,mm0);
% w=diag(1./sum(w,2))*w;
v=det(W2*W2');
j00=1;
k00=1;
iter1=1;
iter2=1;
XX0=XX;
beta=1;
DET1=zeros(3*r,1);
DET2=zeros(3*r,1);
DET=zeros(3*r,1);
WW2=W1'*W1;
 while iter2<3*r     
      iter2;
      tic
       XX=beta*XX0;
    WW=W2*XX*W2';
    A_P=adj_mat(WW);
    i=mod(iter2,r);
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
                if t>=i 
                    tt=t+1; 
                end
                coef10=coef10+(-1)^(t+i-1)*det(temp2)*W2(tt,:)*XX;                
            end
            coef11=coef11+(-1)^(i+j)*XX*W2(j,:)'*coef10;
        end
        if j==i
            coef21=XX*det(temp1);
        end
        if j>i
            coef30=zeros(1,m);
            for t=1:r-1
                 temp3(:,:)=A_P_temp1(t,i,:,:);
                 tt=t;
                if t>=i 
                    tt=t+1;
                end
                coef30=coef30+(-1)^(t+i)*det(temp3)*W2(tt,:)*XX;                
            end
            coef31=coef31+(-1)^(i+j)*XX*W2(j,:)'*coef30;
        end
        obj_coef=coef11+coef21+coef31;
    end

    Aeq=ones(1,m0);
    beq=1;
    H=(obj_coef+obj_coef')/2-aa*WW2*abs(max(max(obj_coef)));
    f=zeros(m,1);
    A=-X2';
    b=zeros(n0,1);
    
    [y,fval1,exitflag1]=quadprog1(-H,f,A,b,Aeq,beq,[],[],W2(i,:)',optimset('Display','off'));

    obj2(1,iter2)=exitflag1; 
    obj2(2,iter2)=fval1; 
    if exitflag1<0
        beta=beta*1.2;
        iter2=iter2-1;        
    end
    

    if exitflag1>-2
        W2(i,:)=y';
    end
    y;
    W2;
   iter2=iter2+1;
     A_est2=inv(W2);
      S_est2=W2*X2; 
        DET2(iter2)=-log(abs(det(S_est2*S_est2')));





 
      
 W=[W1;W2];
 A_est=[A_est1,A_est2];
 S_est12=W1*X2;
 S_est21=W2*X1;
S_est11=[S_est1,S_est12];
S_est22=[S_est21,S_est2];
  S_est=[S_est11;S_est22];     
   Xest=A_est*S_est;    
err(iter2)=norm(X0-Xest,'fro')/sqrt(mm0*nn0);   
  [ratio(iter2)] = dictdist(W,Omega,0.01);
  [recoveryH(iter2)]=dictdist(S_est',H0',0.01);
   DET(iter2)=DET1(iter1)+DET2(iter2);
   toc
end
 
 

beta

%figure
%subplot(2,1,1)
%plot(obj(1,:))
%subplot(2,1,2)
%plot(obj(2,:))
%title('det')

figure;   
plot(DET,'-k');
%title('DET');
   xlabel('Iteration Number');
    ylabel('The sparsity of H'); 

end

function [A S obj]=nLCA_IVM2(X0)
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
index_f=index_s(find(p>0)); % list the indices of data corresponding to possible constraints
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

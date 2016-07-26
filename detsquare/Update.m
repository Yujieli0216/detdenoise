function[Aest Sest Sst Xest Cost Err St ratio]=Update(X,A0,S0,A)
[a,b]=size(S0);
M=sqrt(a);
S=reshape(S0,M,M);
r=M;
Sest=S;
m=size(X,1);
Subm=zeros(m,1);
Stt=zeros(M,1);
Sst=zeros(a,b);
Cost=zeros(M,1);
Err=zeros(M,1);
St=zeros(M,1);
ratio=zeros(M,1);


    SS=S'*S;
    A_P=adj_mat(SS);
   for i=1:M
       tic
    coef11=zeros(M,M);
    coef21=zeros(M,M);
    coef31=zeros(M,M);
    for j=1:M
        temp1(:,:)=A_P(i,j,:,:);
        A_P_temp1=adj_mat(temp1);

        if j<i
            coef10=zeros(1,M);
            for t=1:r-1
                 temp2(:,:)=A_P_temp1(t,i-1,:,:);
                 tt=t;
                if t>=i 
                    tt=t+1; 
                end
                 Stt=S(:,tt);
                coef10=coef10+(-1)^(t+i-1)*det(temp2)*Stt';                
            end
            coef11=coef11+(-1)^(i+j)*S(:,j)*coef10;
        end
        if j==i
            coef21=det(temp1);
        end
        if j>i
            coef30=zeros(1,M);
            for t=1:r-1
                 temp3(:,:)=A_P_temp1(t,i,:,:);
                 tt=t;
                if t>=i;
                    tt=t+1;
                end
                Stt=S(:,tt);
                coef30=coef30+(-1)^(t+i)*det(temp3)*Stt';                
            end
            coef31=coef31+(-1)^(i+j)*S(:,j)*coef30;
        end
    end
        obj_coef=coef11+coef21+coef31;
             p=1+(i-1)*M;
             q=i*M;
             
         for iter=1:2000
             
        Sest(:,i)=Sest(:,i)-0.05*A0(:,p:q)'*(X-Subm)-0.05*0.1*(obj_coef+obj_coef')*S(:,i);
           Si=Sest(:,i);  
            
             
        if min(min(Si)) < 0
           Si(Si < 0) = 0;
        end
       %  if sum(Si) <= 0
       %     Si=Si;
       %  end
         Si=Si./sum(Si);
         Sest(:,i)=Si;   
        A0(:,p:q)=A0(:,p:q)-0.001*(X-Subm)*Sest(:,i)';
         if min(min(A0)) < 0
           A0(A0 < 0) = 0;
         end
        
        Subm=Subsum(A0,Sest);
        Subm=Subm./sum(Subm);
         Err(iter)=sum(sum((X-Subm).^2));
        Cost(iter)=0.5*sum(sum((X-Subm).^2))-0.1*det(Sest'*Sest);
        St(iter)=-0.001*det(Sest'*Sest);
         end
  
      
        
      
        
        Sst=reshape(Sest,a,1);
        Sst=Sst./sum(Sst);
        
         Aest=A0;
         Aest=normcols(Aest);
         Xest=Aest*Sst;
        
        

         [ratio(i)] = dictdist(Aest,A,0.1);
         toc
   end
   
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

    end
end

 



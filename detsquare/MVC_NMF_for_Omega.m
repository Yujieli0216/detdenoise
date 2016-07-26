function [Omega,Xest,St,Err,Cost,ratio] = MVC_NMF_for_Omega(X, Omega, Omega0, H, eps)
[n N] = size(X);
Xest=zeros(size(X));
Err=zeros(1000,1);
St=zeros(1000,1);
Cost=zeros(1000,1);
ratio=zeros(1000,1);
[m M]=size(Omega);
At=[];
Bt=[];
Ct=[];
na=4000;
tic
for iter=1:1000
    tic
   
%if iter>500
%    na=0.05;
%end
%Omega=Omega+0.08*pinv(Omega*X*X'*Omega')*Omega*X*X'-0.002*(Omega*Omega'-eye(m))*Omega;
%Omega=Omega+0.0001*X-0.0001*(Omega*Omega'-eye(m))*Omega;
Omega=Omega-0.00001*(Omega*X-H)*X'+0.00001*na*pinv(Omega*X*X'*Omega')*Omega*X*X';


   if min(min(Omega)) < 1e-6
      Omega(Omega < 1e-6) = 0;

  end
   Omega=normrows(Omega);
  
  Hest=Omega*X;
  Xest=pinv(Omega)*Omega*X;

  
  if min(min(Xest)) < 1e-6
      Xest(Xest < 1e-6) = 0;
  end
   
  Xest=normcols(Xest);
 
  %Err(iter)=sum(sum((X-Xest).^2));
  Err(iter)=sum(sum((H-Hest).^2))/N;
  
  
  At=Omega*X;
  Bt=At*At';
  Bt=normrows(Bt);
  Ct=det(Bt);
  St(iter)=-na*log(abs(Ct))/N;

 

  
 % Cost(iter)=-log(abs(Ct))/log(10)+0.025*0.5*sum(sum((Omega*Omega'-eye(m)).^2));
Cost(iter)=0.5*sum(sum((H-Hest).^2))/N-na*log(abs(Ct))/N;
[ratio(iter)] = dictdist(Omega,Omega0,0.055);




    if rem(iter,50)==0
        disp('done!')
    end
    toc
    
end
Omega=normrows(Omega);
Xest=normcols(Xest);

    
    
    
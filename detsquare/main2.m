
clc
clear
close all
for b=1:1;
% Set problem dimensions
n=5; 
d=n^2;
%d=n^2; % signal dimension
r=15; % subspace dimension for each signal 
N=2000; % number of training signals


% Generate analysis signals lying in r-dimensional nullspaces

% D=rand(d,2*d);
% D=normcols(D);
% D1=D(:,1:d);
% D2=D(:,d+1:2*d);
% Omega1=pinv(D1);
% Omega2=pinv(D2);
% Omega=[Omega1;Omega2];
% D=D./max(max(D));
% h1=figure;
% DisplayS(D,h1);
%  title('D0')




%D1=rand(d,d);

%D=[D1];

D=GenerateOmegaDIF(n);
D=normcols(D);
%D=D./max(max(D));

Omega1=pinv(D);


Omega=[Omega1];

p=size(D,1); 

[X1, H1] = gererateSyntheticDictionaryAndData(D, N, r, 0);

Omega=normrows(Omega);
%D2=X1*pinv(H2);
%D=[D1,D2];


h1=figure;
DisplayOmega(Omega,h1);
 %title('Omega0')
 
 H0=[H1]; 

X=[X1];

%show the sparse of H0
    S=abs(H0)<1e-6;    %the number of zeros
    figure; clf; 
    h=hist(sum(S),1:1:p); 
    hist(sum(S),1:1:p); 
    hold on; 
    plot([d-r d-r],[0 max(h)]); 
    xlabel('cosparsity ');
    ylabel('# of sample signals'); 
    
% 
%  IMAGE=ones(1,d*n+d+1);
%     for j=1:1:n
%         ROW=ones(n,1);
%         for k=1:1:d
%             pos=randperm(N);
%             pos=pos(1);
%             ROW=[ROW, reshape(X(:,pos),n,n),ones(n,1)];
%         end
%         IMAGE=[IMAGE; ROW; ones(1,d*n+d+1)];
%     end
%     figure,imagesc(IMAGE); colormap(gray(256));
%     axis image; axis off;
%  
%Omega2=pinv(D2);
% Omega2=rand(n^2,n^2);
% Omega=[Omega1;Omega2];
% H2=Omega2*X1;



% Omega=rand(2*n^2,n^2);
% %Omega=GenerateOmegaDIF(n);
% [X0,X]=gererateSyntheticDictionaryAndData(Omega,N,r);
% H=Omega*X;








[A_est Omeganew S_est ratio recoveryH DET]=IVM_QP1(X,d,Omega,H0);
S_est=normcols(S_est);




D2=A_est;
D2=normcols(D2);
D2=D2./max(max(D2));
Hest=S_est;

Omeganew=normrows(Omeganew);
h3=figure;
DisplayOmega(Omeganew,h3);
%title('Omega1');


Xest=D2*S_est;    

%show the sparse of H0
    S=abs(S_est)<1e-6;    %the number of zeros
    figure; clf; 
    h=hist(sum(S),1:1:p); 
    hist(sum(S),1:1:p); 
    hold on; 
    plot([d-r d-r],[0 max(h)]); 
    xlabel('cosparsity');
    ylabel('# of recovered signals'); 
    

%  IMAGE=ones(1,d*n+d+1);
%     for j=1:1:n
%         ROW=ones(n,1);
%         for k=1:1:d
%             pos=randperm(N);
%             pos=pos(1);
%             ROW=[ROW, reshape(X(:,pos),n,n),ones(n,1)];
%         end
%         IMAGE=[IMAGE; ROW; ones(1,d*n+d+1)];
%     end
%     figure,imagesc(IMAGE); colormap(gray(256));
%     axis image; axis off;
%    
figure;   
plot(ratio*100,'-k');
%title('recovery');
   xlabel('Iteration Number');
    ylabel('The Rate of Recovery60 (%)'); 
% figure;
% plot(recoveryH,'-k');
% title('recoveryH');
%     
%  
%  

fname1=['resultserr1_smu1_' num2str(b) 'ratio.mat'];
fname2=['resultserr2_smu1_' num2str(b) 'DET.mat'];
save(fname1,'ratio');
save(fname2,'DET');

end


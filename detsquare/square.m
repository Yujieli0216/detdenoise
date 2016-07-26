clc
clear
close all
n=2; 
d=n^2; % signal dimension
r=2; % subspace dimension for each signal 
N=20000; % number of training signals
L1=0;

% Generate analysis signals lying in r-dimensional nullspaces
Omega0 = GenerateOmegaDIF(n); % Omega_DIF
DisplayOmega(Omega0);
%Omega0=rand(n^2,n^2);
p=size(Omega0,1); % number of atoms
%[X,S,L]=GenerateAnalysisSignals(Omega0,d-r,N);
[X, H] = gererateSyntheticDictionaryAndData(Omega0, N, r);
[A_est S_est obj ratio]=IVM_QP(X,d,Omega0);
S_est=normcols(S_est);

%DisplayOmega(S_est);

Omega2=A_est;
Omega2=normrows(Omega2);

Hest=S_est;


h3=figure;
DisplayOmega(Omega2,h3);
title('Omega2')


Xest=Omega2*S_est;    


   
figure;   
plot(ratio,'-k');
title('recovery')
    
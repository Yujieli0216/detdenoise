clear; close all;

s01=imread('data/face_01.tif');
s02=imread('data/face_02.tif');
s03=imread('data/face_03.tif');
s04=imread('data/face_04.tif');
s05=imread('data/face_05.tif');
s06=imread('data/face_06.tif');
[m01 n01]=size(s01);
k1=1;
k2=1;
s1=s01(1:k1:m01,1:k2:n01);
s2=s02(1:k1:m01,1:k2:n01);
s3=s03(1:k1:m01,1:k2:n01);
s4=s04(1:k1:m01,1:k2:n01);
s5=s05(1:k1:m01,1:k2:n01);
s6=s06(1:k1:m01,1:k2:n01);

[m0 n0]=size(s1);

S0=double([reshape(s1,1,m0*n0);reshape(s2,1,m0*n0);reshape(s3,1,m0*n0);reshape(s4,1,m0*n0);reshape(s5,1,m0*n0);reshape(s6,1,m0*n0)]);
S=[S0(1:2,:);S0(3:6,:)];

[n N]=size(S);

A=rand(n,n);

X=A*S;
[n N]=size(X);


[A_est1 H(1:n,:) obj1 ratio]=IVM_QP(X,n,A);%,1e5);

Recovered_S = A_est1*X;


figure;
subplot(2,3,1); imshow(reshape(S(1,:),350,275),[]);
subplot(2,3,2); imshow(reshape(S(2,:),350,275),[]);
subplot(2,3,3); imshow(reshape(S(3,:),350,275),[]);
subplot(2,3,4); imshow(reshape(S(4,:),350,275),[]);
subplot(2,3,5); imshow(reshape(S(5,:),350,275),[]);
subplot(2,3,6); imshow(reshape(S(6,:),350,275),[]);
title('Original Faces');

figure;
subplot(2,3,1); imshow(reshape(X(1,:),350,275),[]);
subplot(2,3,2); imshow(reshape(X(2,:),350,275),[]);
subplot(2,3,3); imshow(reshape(X(3,:),350,275),[]);
subplot(2,3,4); imshow(reshape(X(4,:),350,275),[]);
subplot(2,3,5); imshow(reshape(X(5,:),350,275),[]);
subplot(2,3,6); imshow(reshape(X(6,:),350,275),[]);
title('Mixed Faces');

figure;
subplot(2,3,1); imshow(reshape(H(1,:),350,275),[]);
subplot(2,3,2); imshow(reshape(H(2,:),350,275),[]);
subplot(2,3,3); imshow(reshape(H(3,:),350,275),[]);
subplot(2,3,4); imshow(reshape(H(4,:),350,275),[]);
subplot(2,3,5); imshow(reshape(H(5,:),350,275),[]);
subplot(2,3,6); imshow(reshape(H(6,:),350,275),[]);
title('Unmixed Faces');

figure;   
plot(ratio,'-k');
title('recovery')
    
 
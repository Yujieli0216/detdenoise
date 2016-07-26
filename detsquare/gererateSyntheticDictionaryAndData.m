function [data, coefs] = gererateSyntheticDictionaryAndData(Omega, N, L, SNRdB)


% randn('state',sum(100*clock));
rand('state',sum(100*clock));

Dictionary = Omega; %randn(dim,K);


[data,coefs] = CreateDataFromDictionarySimple(Dictionary, N, L);

if (SNRdB==0) | (SNRdB == 80) 
    return
else
    noise = rand(size(data));  %randn(size(data));
    actualNoise = calcNoiseFromSNR(SNRdB,data, noise);
    SNR = calcSNR(data, data+actualNoise);
    data =  data + actualNoise*SNR/SNRdB;   
end



% %=============================================
% % Display original dictionary 
% %=============================================
% figure;
% displayDictionaryElementsAsImage(Dictionary, 50, 1,1,20,0);
% title('Original dictionary');
% set(gcf,'units','normalized','position',[0,0,1,0.935]);
% 
% %=============================================
% % Display original signal 
% %=============================================
% figure;
% displayDictionaryElementsAsImage(data, 10, 5,5,5,0);
% title('Original signal');
% set(gcf,'units','normalized','position',[0,0,1,0.935]);

function [X,xOrig] = CreateDataFromDictionarySimple(dictionary, numElements, numCoef)
maxRangeOfCoef = 1;
resolution = 0.0001;
n=6;
N=numElements;
[p,d]=size(dictionary);

xOrig = zeros(size(dictionary,2),numElements);
%vecOfValues = -1*maxRangeOfCoef:resolution:maxRangeOfCoef;
%coefs = randsrc(numCoef,numElements,vecOfValues);
coefs = rand(numCoef,numElements)*maxRangeOfCoef;  %randn
xOrig(1:numCoef,:) = coefs;
for i=1:size(xOrig,2)
    xOrig(:,i) = xOrig(randperm(size(xOrig,1)),i);
end
xOrig=normcols(xOrig);
%dictionaryElementIndices = randsrc(numCoef*numElements,1,[1:size(dictionary,2)])   ; 
%matrixOfIndices = repmat([1:numElements],numCoef,1);
%xOrig(sub2ind(size(xOrig),dictionaryElementIndices,matrixOfIndices(:))) = coefs;
X = dictionary*xOrig;

%   if min(min(X)) < 0
%         X(X < 0) = 0;
% 
%    end

X=normcols(X);
S=abs(xOrig)<1e-6;


% Present several examples

    h0=figure;
    DisplayOmega(dictionary,h0);
%     IMAGE=ones(1,20*n+20+1);
%     for j=1:1:10
%         ROW=ones(n,1);
%         for k=1:1:20
%             pos=randperm(N);
%             pos=pos(1);
%             ROW=[ROW, reshape(X(:,pos),n,n),ones(n,1)];
%         end
%         IMAGE=[IMAGE; ROW; ones(1,20*n+20+1)];
%     end
%     figure,imagesc(IMAGE); colormap(gray(256));
%     axis image; axis off;


% Create an histogram of the co-sparsity levels

%     figure(3); clf; 
%     h=hist(sum(S),1:1:p); 
%     hist(sum(S),1:1:p); 
%     hold on; 
%     plot([d-numCoef d-numCoef],[0 max(h)]); 
%     xlabel('cosparsity l');
%     ylabel('# of signals'); 

function  actualNoise = calcNoiseFromSNR(TargerSNR, signal, randomNoise)
signal = signal(:);
randomNoiseRow = randomNoise(:);
signal_2 = sum(signal.^2);
ActualNoise_2 = signal_2/(10^(TargerSNR/10));
noise_2 = sum(randomNoiseRow.^2);
ratio = ActualNoise_2./noise_2;
actualNoise = randomNoiseRow.*repmat(sqrt(ratio),size(randomNoiseRow,1),1);
actualNoise = reshape(actualNoise,size(randomNoise));

function SNR = calcSNR(origSignal, noisySignal)
errorSignal = origSignal-noisySignal;
signal_2 = sum(origSignal.^2);
noise_2 = sum(errorSignal.^2);

SNRValues = 10*log10(signal_2./noise_2);
SNR = mean(SNRValues);


%function [Dictionary, data] = gererateAnalysisDictionaryAndData(N, dim, K)
%rand('state',sum(100*clock));
%Dictionary = rand(K,dim);
%Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
%[data] = CreateDataFromDictionarySimple(Dictionary, N);

%function [D,xOrig] = CreateDataFromDictionarySimple(dictionary, numElements)
%maxRangeOfCoef = 1;
%resolution = 0.0001;
%xOrig = rand(size(dictionary,2),numElements);
%D = dictionary'*dictionary*xOrig;

function [X,S,L]=GenerateAnalysisSignals(Omega,RankNull,N,ShowFigures)
if nargin<4
    ShowFigures=1;
end
[p,d]=size(Omega);
n=sqrt(d);
L=0;

% Choosing Lambda
X=zeros(d,N);



h=waitbar(0,'Gathering signals ...'); 
for k=1:N
    if rem(k,10)==0
        waitbar(k/N);
    end
    List=randperm(p);
    Q=ComputeOrthoSet(Omega(List,:));
    Q=Q(1:RankNull,:);
    x=randn(d,1);
    X(:,k)=(eye(d)-Q'*Q)*x;
end

    


    if min(min(X)) < 0
        X(X < 0) = 0;

   end
close(h)
X=normcols(X);
S=abs(Omega*X)<1e-6;
L=sum(S)/d;

% Present several examples
if ShowFigures && floor(n)==ceil(n)
    h0=figure;
    DisplayOmega(Omega,h0);
    IMAGE=ones(1,20*n+20+1);
    for j=1:1:10
        ROW=ones(n,1);
        for k=1:1:20
            pos=randperm(N);
            pos=pos(1);
            ROW=[ROW, reshape(X(:,pos),n,n),ones(n,1)];
        end
        IMAGE=[IMAGE; ROW; ones(1,20*n+20+1)];
    end
    figure,imagesc(IMAGE); colormap(gray(256));
    axis image; axis off;
end

% Create an histogram of the co-sparsity levels
if ShowFigures
    figure(3); clf; 
    h=hist(sum(S),1:1:p); 
    hist(sum(S),1:1:p); 
    hold on; 
    plot([RankNull RankNull],[0 max(h)]); 
    xlabel('cosparsity l');
    ylabel('# of signals'); 
end

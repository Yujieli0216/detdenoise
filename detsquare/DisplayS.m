function [OmegaDisp] = DisplayS(Omega,h)
[p,d]=size(Omega);
n=round(sqrt(p));
a=min(Omega(:));
b=max(Omega(:));
N1=n;
N2=ceil(d/n);
count=1; 
if nargin<2
    h=figure;
end
  %IMAGE=ones(1,5*n+5+1);

IMAGE=ones(1,1+(n+1)*N2); 
    for j=1:1:N1
        ROW=ones(n,1);
        for k=1:1:N2
            
            ROW=[ROW, reshape(Omega(:,count),[n,n]),ones(n,1)];
            count=count+1;
        end
        IMAGE=[IMAGE; ROW; ones(1,1+(n+1)*N2)];
    end
    
    figure(h);image(IMAGE*255); colormap(gray(256));
    axis image; axis off;
 
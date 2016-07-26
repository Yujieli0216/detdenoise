function [D]=dictionary(m,n)
D=rand([m n]);
D=normcols(D);
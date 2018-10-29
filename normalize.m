function f=normalize(A)
% mean of A for each variable(each column)
n=size(A,1);
mu=mean(A,1); 
% A=A-mean(mu);
% mu=mean(A,1); 
A=A-repmat(mu,n,1);
% variance of each column of A
%  nu=(sum(A.^2,1)/(n)).^(0.5);    
 nu=(n*var(A)).^0.5;
 f=A./repmat(nu,n,1);

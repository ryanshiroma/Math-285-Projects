function [V,dist,W,D,sigma] = nCut(X,tuning_method,tuning_param)
% NCUT  runs the ncut algorithm on X
%   V = NCUT(X) median of the 7th nearest neighbor as sigma.
%   V = NCUT(X,'median',k) median of the kth nearest neighbors as sigma.
%   V = NCUT(X,'average',k) average of the 7th nearest neighbors as sigma.
%   V = NCUT(X,'custom',val) val as sigma.
%   V = NCUT(X,'self-tuning',k) kth nearest neighbor of each node as sigma.

knn=7; %default nearest neighbor setting
n=size(X,1);
if nargin ==1
    tuning_method = 'median';
elseif strcmp(tuning_method,'custom')
    sigma=ones(n,1)*tuning_param;
end

if ~strcmp(tuning_method,'custom')
%Lets first initialize sigma
    dist=zeros(n,n); % initialize the distance matrix
    for i=1:n
        for j=i+1:n
            dist(i,j)=sqrt(sum((X(i,:)-X(j,:)).^2));
        end
    end
    dist=dist + dist'; 
    sorted_dist=sort(dist,2);
    
    if strcmp(tuning_method,'average') 
        sigma=ones(n,1)*mean(sorted_dist(:,knn+1));
        
    elseif strcmp(tuning_method,'median') 
        sigma=ones(n,1)*median(sorted_dist(:,knn+1));
         
    elseif strcmp(tuning_method,'self-tuning')
        sigma=sorted_dist(:,knn+1);
    end
end 
sigma
%now lets creat the W matrix
W=zeros(n,n);
for i=1:n
    % Since K is symmetric, we only need to compute an upper triangular
    % matrix just add the transpose to itself.
    for j=i:n 
        W(i,j)= exp( -sum((X(i,:)-X(j,:)).^2)/(2*sigma(i)*sigma(j)));
    end
end
W=W + W'-2*diag(diag(W));% compute W from the upper triangular W matrix.

D=diag(sum(W)); %creates the D matrix
L_rw=eye(n)-inv(D)*W; %creates the L random walk matrix
[V,eigenvalues]=eig(L_rw); %finds the eigenvectors/values

[~,index]=sort(diag(eigenvalues)); %sorts the eigenvectors by eigenvalue.
V=V(:,index);
return;
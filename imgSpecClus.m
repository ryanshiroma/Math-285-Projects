function Y = imgSpecClus(Ximg,clusters,sigma)
   Ximg= imrez(Ximg,50,50);
X = zeros(50*50,3);
double(Ximg);
count = 1;
for i = 1:50
    for j = 1:50
        X(count,1) = i;
        X(count,2) = j;
        X(count,3) = Ximg(i,j,1);
        count = count +1;
    end
end 

N = size(X,1);

    avgDist = sigma;
%Now the weighted matrix
W = zeros(N);
for i = 1:N
    for j = 1:N
        if i ~= j 
            if (abs(X(i,1) - X(j,1))) <10 || (abs(X(i,2) - X(j,2))) <10 
                W(i,j) = exp(-(norm(X(i,:)-X(j,:))^2)/(2*avgDist^2));
            end
        end
    end
end
%then the D matrix
D = diag(sum(W,2));
%Then either the random walk matrix or unormalized
   % L = eye(N) - inv(D)*W;
    L = D-W;
[U,S] =eig(L);
[S,I] = (sort(diag(S),'ascend'));
U= U(:,I);
[Y,~,~]= kmeans(U(:,2:clusters),clusters,'Replicates',10);
end
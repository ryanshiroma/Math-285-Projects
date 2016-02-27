function A= superPixelAffinity(X,l,Am,C)

spixels=size(unique(l),1);
binranges = -10:30:265;
A=zeros(spixels,spixels); % initialize the distance matrix
distr=zeros(spixels,size(binranges,2));
for i=1:spixels
    hist_c_out=[histc(X(l==i),binranges)];
    if size(hist_c_out,2)>1
        hist_c_out=hist_c_out'
    end
    distr(i,:)= hist_c_out./[repmat(sum(sum(l==i)),size(binranges,2),1)]; 
end

for i=1:spixels
    for j=i+1:spixels
        if Am(i,j)==1
            BC=-log(sum(sqrt(distr(i,:).*distr(j,:))));
            %me
            %HD=-1/sqrt(2)*(sum(sqrt(distr(i,:))-sqrt(distr(j,:)))^2);
 
            A(i,j)=exp(-BC);
           %intensity=0.02*(mean(X(l==i))-mean(X(l==j)))^2;
           %intensity
           %c=[C.c];
           %r=[C.r];
           %d=0.005*((r(i)-r(j))^2+(c(i)-c(j))^2);
           %d
           %A(i,j)=exp(-d-intensity);
        else
            A(i,j)=0;
        end
    end
end
size(A)
'done'
A=A+A';
end
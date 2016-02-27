function A= imageAffinity(X)

pixels=size(X,1)*size(X,2);
n1=size(X,2);
n2=size(X,1);
A=zeros(pixels,pixels); % initialize the distance matrix
for i=1:pixels
    x1=floor(((i-1)/n1)+1);
    y1=mod(i-1,n2)+1;
    for j=i+1:pixels
        x2=floor(((j-1)/n1)+1);
        y2=mod(j-1,n2)+1;
        if abs(x2-x1)<=5 && abs(y2-y1)<=5;
        %if sqrt((x2-x1)^2+(y2-y1)^2)<=4;
            A(i,j)=exp(-(((X(x1,y1)-X(x2,y2))^2)*.005+((x1-x2)^2+(y1-y2)^2)*0.6));
        end
    end
end
A=A+A';
end
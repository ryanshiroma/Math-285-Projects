function Y = imageKmeans(Ximg,k)
    %Applys kmeans to an image Ximg
    %before running kmeans it restructers the 
    %image so that x and y positions are in the first
    %two columns and red intensities are in third
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
    Y = kmeans(X,k,'Replicates',10);
end
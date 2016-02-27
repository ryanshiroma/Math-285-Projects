%%% SCRIPT FILE FOR THE MATH 285 PROJECT



%%%%%%%%%%%%% Let's first test the clustering on some simple data sets %%%%%%%%

%%%% Spectral Clustering on the circle images
circles=rgb2gray(imread('circles.png'));
circles(circles>50)=255;
circles(circles<=50)=10;
circles=circles(1:15:600,1:15:600);
imagesc(circles);
n=size(circles,1);
X=reshape(circles,n*n,1);
X=double(circles);

%%%% Spectral Clustering on the circle images
shapes=rgb2gray(imread('shapes.png'));
shapes=shapes(1:15:300,1:15:400);
shapes=shapes(1:20,1:20);
n=size(shapes,1);
X=reshape(shapes,size(shapes,1)*size(shapes,2),1);
X=double(shapes);



%%%%%%% Now let's run our spectral clustering on the real video %%%%%%%
v = VideoReader('IMG_3304.mp4');
vid=zeros(192,341,size(1:6:400,1));
for i=1:6:400
 temp=read(v,i);
vid(:,:,ceil(i/6))=temp(:,:,1);
end

clusters=2; %% specify only two clusters for now (hand and background)
vid=read(v,1); %%% test on the first frame of the video
vid2=vid(:,150-80:end-80,1);
X2=vid2(:,:,1);
X=X2(1:uint8(1/3*end),1:uint8(1/3*end));
imagesc(X);axis equal;
X=double(X);

%%%% Whats the run time on this single frame? %%%%%
tic
A=imageAffinity(X);
V=SpectralClustering(A,clusters);
toc

%%% the runtime is too long, lets subsample every 30th frame %%%%
X2=vid(:,:,1);
imagesc(X2);
X=X2(1:30:end,1:30:end);
X=double(X);
n1=size(X,1);
n2=size(X,2);
imagesc(X)

%%% every 10th frame %%%%%
X=X(1:10:end,1:10:end,1);
n1=size(X,1);
n2=size(X,2);
Xnew=zeros(n1*n2,3);

count=1
for i=1:n1
    for j=1:n2
        Xnew(count,1)=i;
        Xnew(count,2)=j;
        Xnew(count,3)=X(i,j);
        count=count+1;
    end
end

plot(Xnew,'.')
labels_kmeans=kmeans(Xnew,3,'Replicates', 10);
labels=reshape(labels_kmeans,n2,n1)';
imagesc(labels);


%%%% Let's cycle through the movie file frame by frame %%%%
for i=1:size(vid,3)      

    pause(0.1);
    imagesc(vid(:,:,i));

 
end

imagesc(shapes)
pixels=size(X,1)*size(X,2);




%%%%% Let's look for a better solution. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SLIC Superpixels reimplementation
%% http://infoscience.epfl.ch/record/149300



load Hand;

clusters=2;

[l, Am, C] = slic(temp, 30, 15, 1, 'median');

imagesc(l)
X(:,:,2)=100*ones(1080,1080);
X(:,:,3)=100*ones(1080,1080);
imshow(drawregionboundaries(l, X, [255 255 255]))
    




v = VideoReader('IMG_3304.mp4');
vid=zeros(192,341,size(1:6:400,1));
for i=1:6:400
 temp=read(v,i);
vid(:,:,ceil(i/6))=temp(:,:,1);
end

vid2=vid(:,150-80:end-80,:);
X=vid2(:,:,1);
X=X(1:10:end,1:10:end);
imagesc(X);axis equal;


width=size(vid2(:,:,1),1);
f_seg=zeros(n,n);
  base_dist1= zeros(255,1);
  base_dist2= zeros(255,1);
        
for p=1:50
tic
n=16
X=vid2(:,:,p);
X_reduced=X(1:10:end,1:10:end);
A=imageAffinity(X);

_reduced);

V=SpectralClustering(A,clusters);

n1=size(X_reduced,1);% 8
n2=size(X_reduced,2);% 8

labels_kmeans = kmeans(V(:,2:clusters), clusters, 'Replicates', 10);
labels=reshape(labels_kmeans,n1,n2)';

        
base_dist1= histc(X_reduced(labels==1),1:255)./repmat(sum(sum(labels==1)),255,1);
base_dist2= histc(X_reduced(labels==2),1:255)./repmat(sum(sum(labels==2)),255,1);

for i=1:(width/n)
    for j=1:(width/n)
     
        X_block=X(n*(i-1)+1:n*i,n*(j-1)+1:n*j);
        imagesc(X_block)
        A=imageAffinity(X_block);

        V=SpectralClustering(A,clusters);

        n1=size(X_block,1);% 8
        n2=size(X_block,2);% 8

        labels_kmeans = kmeans(V(:,2:clusters), clusters, 'Replicates', 10);
        labels=reshape(labels_kmeans,n1,n2)';

        %initialize base_dist


        % set up new distributions
        new_dist1= histc(X_block(labels==1),1:255)./repmat(sum(sum(labels==1)),255,1);
        new_dist2= histc(X_block(labels==2),1:255)./repmat(sum(sum(labels==2)),255,1);

        %assign cluster to 'label_1' cluster

        BC_1base1 =sum(sqrt(base_dist1.*new_dist1));
        BC_1base2 =sum(sqrt(base_dist2.*new_dist1));
        BC_2base1 =sum(sqrt(base_dist1.*new_dist2));
        BC_2base2 =sum(sqrt(base_dist2.*new_dist2));
        
        
        %four possibilities
        if BC_1base1>BC_1base2 && BC_2base2>BC_2base1 %%normal
             f_seg(n*(i-1)+1:n*i,n*(j-1)+1:n*j)=(labels==1);
             base_dist1= ((i*j)*base_dist1+new_dist1)/(i*j+1);
             base_dist2= ((i*j)*base_dist2+new_dist2)/(i*j+1);
        elseif BC_1base1<BC_1base2 && BC_2base2<BC_2base1 %%reversed
             f_seg(n*(i-1)+1:n*i,n*(j-1)+1:n*j)=(labels==2);
             base_dist1= ((i*j)*base_dist1+new_dist2)/(i*j+1);
             base_dist2= ((i*j)*base_dist2+new_dist1)/(i*j+1);
        elseif BC_1base1>BC_1base2 && BC_2base1>BC_2base2 %% 1 cluster
             f_seg(n*(i-1)+1:n*i,n*(j-1)+1:n*j)=ones(n,n);
             base_dist1= ((i*j)*base_dist1+(new_dist2+new_dist1)/2)/(i*j+1);
        elseif BC_1base1<BC_1base2 && BC_2base1<BC_2base2  %% 1 cluster reverse
            f_seg(n*(i-1)+1:n*i,n*(j-1)+1:n*j)=zeros(n,n);
             base_dist2= ((i*j)*base_dist2+(new_dist2+new_dist1)/2)/(i*j+1);
        end
          
    end
end
toc
imagesc(f_seg)
drawnow;
end





%%%%% movie
for i=1:size(vid,3)      

    pause(0.1);
    imagesc(vid(:,:,i));
end



vid=zeros(192,341,3,size(1:6:400,1));
%%%extract color
for i=1:6:400
    vid(:,:,:,i)=read(v,i);
end


clear all;
clear v;
v = VideoReader('IMG_3304-2.mp4');

clusters=2;
superpixels=2;
for f=1:50
    
    X=read(v,f*2);
    X=uint8(X);
    tic
    [l, Am, C] = slic(X, superpixels, 30, 1, 'median');  
toc
    A=superPixelAffinity(X(:,:,1),l,Am,C);

    V=SpectralClustering(A,clusters);
    n1=size(X,1);% 8
    n2=size(X,2);% 8
    labels_kmeans = kmeans(V(:,2:clusters), clusters, 'Replicates', 10);
    labels=zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            labels(i,j)=labels_kmeans(l(i,j));
        end
    end
    labels(labels==2)=X(labels==2);
    imagesc(labels);
    drawnow;
end





clear all;
clear v;
v = VideoReader('IMG_3304.mp4');
vid=zeros(192,341,:,size(1:6:400,1));



v = VideoReader('IMG_3304-2.mp4');

X=vid2(:,:,:,200);
clusters=2;
superpixels=20;

    X=read(v,265);
  
    X=uint8(X);
    imagesc(X(:,:,1))
    tic
    [l, Am, C] = slic(X, superpixels, 30, 1, 'median');  
    toc
    n1=size(X,1);% 8
    n2=size(X,2);% 8
    labels=zeros(n1,n2);
    avg_int=zeros(55,1);
    for i=1:55
        avg_int(i)=mean(nonzeros(mean(nonzeros(X(l==i)))))
    end
    for i=1:n1
        for j=1:n2
            labels(i,j)=avg_int(l(i,j));
        end
    end
    imagesc(labels)
    imagesc(l)
    
               c=[C.c]
       r=[C.r]  
   
 %plot(c,r,'o');
  
       imagesc(X(:,:,1)); hold on;
       %colormap([1 1 1; 0 0 0]); hold on; 
        imagesc(edge(l,'Canny',0));hold on;
   for i=1:length(C)
   
       for j=1:length(C)
           if Am(i,j)==1
               'test'
              plot([c(i) c(j)],[r(i) r(j)],'go', 'Markersize',10);hold on;
           end
       end
   end
   
   imshow(drawregionboundaries(l, X, [255 255 255]))
  
   
   imagesc(X);hold on;
   for i=1:7
       for j=1:7
           plot(i,j,'go', 'Markersize',10);hold on;
       end
   end
   
   
    %if sum(labels==1)>sum(labels==2)
    %    labels=-labels+repmat(3,size(labels,1),1);
    %end
    labels(labels==2)=X(labels==2);

   
           c=[C.c]
       r=[C.r]  
  % colormap([1 1 1; 0 0 0.5]); hold on; 
       imagesc(labels);hold on;
  % imagesc(edge(l,'Canny',0));hold on;
   for i=1:length(C)
       for j=1:length(C)
           if Am(i,j)==1
               'test'
              plot([c(i) c(j)],[r(i) r(j)],'wo-','LineWidth',2, 'Markersize',15);hold on;
           end
       end
   end
   
   imagesc(X(:,:,1))
   
plot([c(i) c(j)],[r(i) r(j)],'o-')

    imagesc(Q)
    imagesc(X(:,:,1));
    imagesc(X(1:7:end,1:7:end,1))
    drawnow;

    X=(X(1:156:end,1:156:end,1))
    X=vid2(:,:,:,200);
clusters=2;
superpixels=200;
    
    tic
    [l, Am, C] = slic(X, superpixels, 5, 1, 'median');  
toc
tic
    A=superPixelAffinity(X(:,:,1),l,Am,C);
    toc
A=imageAffinity(double(X));
tic
    V=SpectralClustering(A,clusters);
    toc
    toc
    n1=size(X,1);% 8
    n2=size(X,2);% 8
    
    labels_kmeans = kmeans(V(:,2:clusters), clusters, 'Replicates', 10);
     labels=reshape(labels_kmeans,n1,n2)';
    colormap([1 1 0; 0 0 0.5]);
    imagesc(labels);
        
    tic
    kmeans(X3(1:1000),2, 'Replicates', 10);
    toc
    labels=zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            labels(i,j)=labels_kmeans(l(i,j));
        end
    end

    colormap([1 1 0; 0 0 0.5]);
    imagesc(labels);
    imagesc(X)
    drawnow;

    A2=A([find(labels_kmeans==1)' find(labels_kmeans==2)'],:)
    
    

    A3=A2(:,[find(labels_kmeans==1)' find(labels_kmeans==2)'])
    
    
    
    binranges = 1:1:265;
    X(l==2)./[repmat(50390,27,1)]
    for i=1:15
        i=9
    hist(double(X(l==i)),binranges)
    xlim([0 255]);
    pause(1)
    end
    imagesc(l)
    ./[repmat(sum(sum(l==2)),size(binranges,2),1)],binranges))
    hist_c_out=[histc(X(l==2),binranges)];
    
    
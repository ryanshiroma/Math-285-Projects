function [tdx tdy tdz] = projection3D(init_width,x,y,cluster_pix,x_dim,y_dim)

%init_width   is the width of the object in real units(inch, cm, etc)
%x            is the image x location
%y            is the image y location
%cluster_pix  is the total number of pixels in the cluster
%x_dim        is the width in pixels of the image
%y_dim        is the height in pixels of the image


% equivalent focal length on the galaxy phone is 16mm which corresponds to
% a 96.7 degree field of view

%%%% 96.7 = 2*atand(36/(2*16)) equation for this calculation

tot_pixels=x_dim*y_dim; %get total area of the image
%%%%% Z axis
theta = asind(sqrt(cluster_pix/tot_pixels)*sind(96.7/2)); %calculate the angle field of view that the cluster takes up
tdz=init_width/tand(theta);  %calculate the z axis coordinate

%%%x and y axes

theta = asind(((x-(x_dim/2))/x_dim/2)*sind(96.7/2)); %calculate the angle field of view that the cluster takes up
tdx=tand(theta)*tdz;  %calculate the z axis coordinate

theta = asind(((y-(y_dim/2))/y_dim/2)*sind(96.7/2)); %calculate the angle field of view that the cluster takes up
tdy=tand(theta)*tdz;

% tdx=tand((96.7/2)*(x-(x_dim/2)))*tdz; %tan(degrees to the left or right) * z axis
% tdy=tand((96.7/2)*(y-(y_dim/2)))*tdz; %tan(degrees to the up or down) * z axis
end
clear
close all
clc
%%%%%%%% Parameters for reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1_max   = 240;                     % Parameter for reconstruction Matrix
r2_max   = 240;
r3_max   = 240;
N=2*r1_max+1;
N1        = 2*r1_max+1;             % 3D image side length
N2        = 2*r2_max+1;
N3        = 2*r3_max+1;    
dims     = [N1,N2,N3];
u_max    = 100;                     % Projection side length
nuv      = [1,1];                   % Number of subpixels, see buildSystemMatrix
vpRatio  = 0.5;                     % Voxel-to-pixel ratio, see buildSystemMatrix
num_proj = 36;                      % Number of projections to reconstruct from
angle_phiti_interval=360/num_proj;
rnl      = 0.01;                    % Relative noise level
thetaZ=30;                          % Angle for the projection with z

%%%%%% Choose a number of random projection directions  %%%%%%%%%%%%%%%%

listZ=cos(thetaZ/180*pi);
listR_xy=sin(thetaZ/180*pi);
for ii=1:num_proj
    v_list(ii,:)=[listR_xy*cos(2*pi/num_proj*(ii-1)),listR_xy*sin(2*pi/num_proj*(ii-1)),listZ];
end        

%%%%%% Set up the parallel beam system matrix  %%%%%%%%%%%%%%%%%%%%%%%%%

[A,p_all] = buildSystemMatrix(r1_max,u_max,v_list,nuv,vpRatio);


%%%%%% pretreat the emperimental data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_project=zeros(2*u_max+1,2*u_max+1,num_proj);

path='...\';
cd(path);
Sample = matfile('Proj_data.mat');
project = Sample.distort_data;
for n=1:num_proj
v1=v_list(n,:);
phit   = atan2(v1(2),v1(1));
Im=imrotate(project(:,:,n),180-phit*180/pi,'crop');
im_project(:,:,n)=Im;
end
b=im_project(:);

%% Reconstruct
tol = 1e-6;
maxit = 30;
x_sol = lsqr(A,b,tol,maxit);
sig_max=max(x_sol);
s=(find(x_sol<sig_max/10)) ;
x_sol(s)=0;
X_sol = reshape(x_sol,dims);

for i=1:size(X_sol,3)
    X_sol(:,:,i)=1.5*(1-i/size(X_sol,3))^2*X_sol(:,:,i);
end 

X_sol=X_sol/max(max(max(X_sol)));

%%%%%% Display  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N3
    figure(11); imagesc(max(X_sol(:,:,i),0)); axis off; axis equal; colormap gray; colorbar;
    title(num2str(i));
    pause(0.1);
end


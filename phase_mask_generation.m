clear all; close all

%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%
        
        width  =  1920;   % width of original mask
        height =  1080;   % height of original mask
        H = 500;          % height of the side-lobe suppresed mask
               
%%%%%%%%%%%% generating phase mask %%%%%%   
%
pixel_y0 = 0;
z = 1e+6;
for i = 1:1080
    for j = 1:1080       
        pixel_x(i,j) = (-540 + i);
        pixel_y(i,j) = (-540 + j);
        new_pixel_x(i,j)=cos(45/180*pi)*pixel_x(i,j)+sin(45/180*pi)*pixel_y(i,j);
        new_pixel_y(i,j)=cos(45/180*pi)*pixel_y(i,j)-sin(45/180*pi)*pixel_x(i,j);
        pixel_val(i,j) = 0.3*((new_pixel_x(i,j))^3 + (new_pixel_y(i,j))^3) ...
            - 0*(1*(pixel_x(i,j))^2 + 1*(pixel_y(i,j)-pixel_y0)^2) ...
            - 0 * pixel_y(i,j)^2;          
        pixel_val(i,j) = pixel_val(i,j)/z;          
        if (abs(pixel_val(i,j)) > (2*pi))
            pixel_val(i,j) = mod(pixel_val(i,j), (2*pi));    
        else if (pixel_val(i,j) < 0)
            pixel_val(i,j) = (2*pi) + pixel_val(i,j); 
            end
        end     
        phasemask(i,j) = pixel_val(i,j)/2/pi;  
    end
end
figure(1);
imshow(phasemask);

%%%%%%%%%%%% side-lobe suppressed mask %%%%%%

phasemask_b=phasemask(291:790,:);
figure(2);
imshow(phasemask_b);

%%%%%%%%%%%% side-lobe suppressed mask with chirp %%%%%%

lambda = 680e-9;
focus = 20e-2;
fsf = 8e-6;
dx=1;
dy=1;
z0 =5e-3;
[chirpX, chirpY] = meshgrid(((-539:1:540)+dx)/(lambda*focus/fsf),((-539:1:540)+dy)/(lambda*focus/fsf));
chirp_angle = sqrt(max(0,lambda^(-2)-chirpX.^2-chirpY.^2));
chirp_mask = phasemask*255/256*2*pi-2*pi*chirp_angle*z0;
angle_total= (chirp_mask*256/2/pi-256*floor(chirp_mask/2/pi))/256;
phasemask_c=angle_total(291:790,:);
figure(3);
imshow(phasemask_c);
%}
%%%%%%%%%%%% side-lobe suppressed mask with chirp with grating %%%%%%

 Grat_period=6;         % period of grating
 for n=1:width          % generating grating
    pixel(n)=mod(width-n,Grat_period)/Grat_period;
end  
Grating=repmat(pixel,height,1); 
phasemask_d=phasemask_c+Grating(291:790,421:1500);
for i=1:500
    for j=1:1080
        if phasemask_d(i,j)>1 phasemask_d(i,j)=phasemask_d(i,j)-1;end
    end
end
figure(4);
imshow(phasemask_d);

%%%%%%%%%%%%%%%%%%%rotation for S2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Grat_period=50;        % period of grating
 for n=1:width          % generating grating
    pixel(n)=mod(width-n,Grat_period)/Grat_period;
end  
Grating=repmat(pixel,height,1); 

theta=45;               %%%% angel for rotation

    for i = 1:1530
    for j = 1:1530
        
        pixel_x(i,j) = (-765 + i);
        pixel_y(i,j) = (-765 + j);

        new_pixel_x(i,j)=cos((theta+45)/180*pi)*pixel_x(i,j)+sin((theta+45)/180*pi)*pixel_y(i,j);
        new_pixel_y(i,j)=cos((theta+45)/180*pi)*pixel_y(i,j)-sin((theta+45)/180*pi)*pixel_x(i,j);
        pixel_val(i,j) = 0.17*((new_pixel_x(i,j))^3 + (new_pixel_y(i,j))^3) ...
            - 0*(1*(pixel_x(i,j))^2 + 1*(pixel_y(i,j)-pixel_y0)^2) ...
            - 0 * pixel_y(i,j)^2;   
       
        pixel_val(i,j) = pixel_val(i,j)/z;
         
    
        if (abs(pixel_val(i,j)) > (2*pi))
            pixel_val(i,j) = mod(pixel_val(i,j), (2*pi));    
        else if (pixel_val(i,j) < 0)
            pixel_val(i,j) = (2*pi) + pixel_val(i,j); 
            end
        end
        
        phasemask(i,j) = pixel_val(i,j)/2/pi;
   
    end
    end
    %%%%%%%%%%%%%%%cut the mask%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phase_1=zeros(1920,1920);
    phase_1(421:1500,421:1500)=phasemask(226:1305,226:1305);
    phase_1=phase_1(421:1500,:);

for i=1:height
    for j=1:width
        if phase_1(i,j)>1 phase_1(i,j)=phase_1(i,j)-1;
        end
    end
end

 pixel_val2=ones(500,1530);
 piexl_val3= imrotate(pixel_val2,theta);
 [sizex,sizey]=size(piexl_val3);
 if sizex>1080
     x1=floor((sizex-1080)/2);
    
     piexl_val3=piexl_val3(x1:x1+1080-1,:);
 end
 if sizey>1080
     y1=floor((sizey-1080)/2);
     piexl_val3=piexl_val3(:,y1:y1+1080-1);
 end
 
phasemask=piexl_val3;
image_grating1 = zeros(1080,1920);
[maskx,masky]=size(phasemask);
    XX=floor((1080-maskx)/2)+1;
    YY=floor((1920-masky)/2)+1;
image_grating2=image_grating1(XX:XX+maskx-1,YY:YY+masky-1);
SS=image_grating2.*(~logical(phasemask))+phasemask;
image_grating1(XX:XX+maskx-1,YY:YY+masky-1)=SS;

image_pm = zeros(1080,1920);
image_pm =image_grating1;
image_pm_not=not(image_pm);
image_2=image_pm_not.*Grating;

phasemask=phase_1.*image_pm+image_2;

figure(5);
imshow(phasemask);


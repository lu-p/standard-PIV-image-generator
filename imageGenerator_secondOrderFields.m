clear
close all
clc

% output for each image pair: Ifin1, Ifin2, Vfin 

mkdir inputData
mkdir targetData

n_imagepairs=2; % number of image pairs

extrapx=4; % number of pixels beyond every boundary (particles are allowed to leave/enter the camera field of view)

% resx=resy must be set (r0)
resx=40; % image resolution [pixels]
resy=32;

xmax=2*extrapx+resx;
ymax=2*extrapx+resy;

xvf=-xmax/2+0.5:1:xmax/2-0.5; % pixels along x (x axis: sx->dx)
yvf=ymax/2-0.5:-1:-ymax/2+0.5;% pixels along y (y axis: down->top)

C=zeros(resy,resx,2);

for nip=1:n_imagepairs % for each image pair

I1=zeros(ymax, xmax); % first frame
I2=zeros(ymax, xmax); % second frame

% LOCATION OF PARTICLE CENTRES
rho=0.029; %[ppp]
n_particles=round(rho*xmax*ymax); % number of particles to generate frames

r0=[xmax/2*(-1 + 2.*rand(n_particles,1)), ymax/2*(-1 + 2.*rand(n_particles,1))];% initial particle centres (x0, y0)[pixel], (xmax=2*extrapx+resx)

% 2D VELOCITY FIELD
u0=(-4+8*rand); % [px/frame] from -4 to 4
v0=(-4+8*rand);
%u0=0; % uncomment to semplify velocity fields
%v0=0; % uncomment to semplify velocity fields

Ju=(-0.05+0.1*rand(1,2)); % [1/frame] from -0.05 to 0.05 
Jv=(-0.05+0.1*rand(1,2));
%Ju=zeros(1,2); % uncomment to semplify velocity fields
%Jv=zeros(1,2); % uncomment to semplify velocity fields

%Hu=(-0.001+0.002*rand(2,2));%(-0.001 + 0.002*rand(2,2)); % [1/(px*frame)] from -0.001 to - 0.001
%Hv=(-0.001+0.002*rand(2,2));%(-0.001 + 0.002*rand(2,2)); % [1/(px*frame)]
Hu=zeros(2,2); % uncomment to semplify velocity fields
Hv=zeros(2,2); % uncomment to semplify velocity fields


dec=4; % number of significant decimals
u0=round(u0,dec);
v0=round(v0,dec);
Ju=round(Ju,dec);
Jv=round(Jv,dec);
Hu=round(Hu,dec);
Hv=round(Hv,dec);


uu= @(xx,yy) u0+Ju*[xx;yy]+[xx,yy]*Hu*[xx;yy]; % mathematical structure of velocity fields
vv= @(xx,yy) v0+Jv*[xx;yy]+[xx,yy]*Hv*[xx;yy];

Vfint=[u0, v0, Ju, Jv Hu(1,1) Hu(1,2) Hu(2,1) Hu(2,2) Hv(1,1) Hv(1,2) Hv(2,1) Hv(2,2)]; % complete velocity field
Vfin=Vfint(1:6); % first order approximation of velocity fields

s=sprintf('target_%d',nip);
save(fullfile('targetData',s),'Vfin') % first approximation of velocity field is saved as target


% INTEGRATION OVER TIME (multiplication by time interval since fields are stationary)

h=0.5; % half time step


I0=randi(56,[n_particles,1])+199; % peak intensity % [grey value] 200 - 255

dp=1.5+2.5*rand(n_particles,1); % effective particle image diameter [px] from 1.5 to 3.5

for np=1:n_particles
    
    % particle centres, first exposition
    rout1u=r0(np,1)-h*uu(r0(np,1),r0(np,2)); %horizontal component
    rout1v=r0(np,2)-h*vv(r0(np,1),r0(np,2)); %vertical component

    % particle centres, second exposition
    rout2u=r0(np,1)+h*uu(r0(np,1),r0(np,2)); %horizontal component
    rout2v=r0(np,2)+h*vv(r0(np,1),r0(np,2)); %vertical component

for i=1:length(yvf)
     for j=1:length(xvf)
 I1(i,j)=I1(i,j)+I0(np)*exp((-(xvf(j)-rout1u)^2-(yvf(i)-rout1v)^2)/(1/8*dp(np)^2));
 I2(i,j)=I2(i,j)+I0(np)*exp((-(xvf(j)-rout2u)^2-(yvf(i)-rout2v)^2)/(1/8*dp(np)^2));
     end
end
end

figure % see particle images before adding noise and cropping
imshow(I1, 'initialmagnification',1000)
figure
imshow(I2,'initialmagnification',1000)

% crop images to allow particles to leave and enter the field of view
I1=imcrop(I1,[extrapx+1,extrapx+1,resx-1,resy-1]); % [XMIN YMIN WIDTH HEIGHT]
I2=imcrop(I2,[extrapx+1,extrapx+1,resx-1,resy-1]); % [XMIN YMIN WIDTH HEIGHT]

I1(I1>255)=255; % best practice
I2(I2>255)=255;

div_f=255; % div=255: images are normalized in [0,1], div=1: images are in [0, 255]
I1=round(I1/div_f,dec);
I2=round(I2/div_f,dec);

noise1=imnoise(I1,'gaussian',0,max(max(I1))/100); % add noise
noise2=imnoise(I2,'gaussian',0,max(max(I2))/100);

I1=I1+noise1;
I2=I2+noise2;


C(:,:,1)=I1;
C(:,:,2)=I2;
s=sprintf('input_%d',nip);
save(fullfile('inputData',s),'C') % save particle images

figure
imshow(I1, 'initialmagnification',1000)
figure
imshow(I2,'initialmagnification',1000)



save dataset_features
end
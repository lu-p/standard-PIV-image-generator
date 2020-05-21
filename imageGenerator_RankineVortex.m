clear
close all
clc


% output: Ifin1, Ifin2, ufin, vfin

% this script simulates PIV experiments on stationary velocity fields:
% the same velocity field is used to generate all particle image pairs
% Particle image pairs differ for the position of particles (random)

n_imagepairs=2; % number of particle image pairs
dec=4; % significant decimals
extrapx=4; % number of pixels beyond each boundary (particles are allowed to leave/enter the camera field of view)

resx=400; % image resolution x [px]
resy=450; % image resolution along y [px]
xmax=2*extrapx+resx; % image resolution along x before cropping
ymax=2*extrapx+resy; % image resolution along y before cropping

xvf=-xmax/2+0.5:1:xmax/2-0.5; % pixels along x (x axis: sx->dx)
yvf=ymax/2-0.5:-1:-ymax/2+0.5;% pixels along y (y axis: down->top)
[XVF, YVF]=meshgrid(xvf, yvf);

rho=0.029; %[ppp] 
n_particles=round(rho*xmax*ymax);% number of particles before cropping

Ifin1=zeros(n_imagepairs,resy,resx); % first exposition
Ifin2=zeros(n_imagepairs,resy,resx); % second exposition

%% 2D VELOCITY FIELD: RANKINE VORTEX

% vortex (central vortex)
v0=4;%-4+8*rand; %max velocity
[x,y]=meshgrid(xvf,yvf);
[o,r]=cart2pol(x,y);
R1=300;%min(xmax/2,ymax/2)*rand; %radius
uoin = (r <= R1);
uout = (r > R1);
uo = uoin+uout;
uo(uoin) =  v0*r(uoin)/R1;
uo(uout) =  v0*R1./r(uout);
uo(isnan(uo))=0;
u = -uo.*sin(o);
v = uo.*cos(o);
ufin=imcrop(u,[extrapx+1,extrapx+1,resx-1,resy-1]);
vfin=imcrop(v,[extrapx+1,extrapx+1,resx-1,resy-1]);

h=0.5; % half time step

div_f=255; % div_f=255 normalizes images in [0, 1] (real numbers)

I0=randi(56,[n_particles,n_imagepairs])+199; % peak intensity % [grey value] from 200 to 255
dp=1.5+2.5*rand(n_particles,n_imagepairs);% effective particle image diameter [px]

str_imp=num2str(n_imagepairs);
name=sprintf('ig_rkv_%dx%d_%s_noise_rhocost_div%d_dec%d_R1=%d_ok',resx,resy,str_imp, div_f, dec, R1)


for nip=1:n_imagepairs
tic
nip
I1=zeros(ymax,xmax);
I2=zeros(ymax,xmax);

% LOCATION OF PARTICLE CENTRES

r0=[xmax/2*(-1 + 2*rand(n_particles,1)), ymax/2*(-1+2*rand(n_particles,1))];% initial particle centres

% check true velocity field
figure
quiver(x,y,u,v)
xlabel('x [px]')
ylabel('y [px]')
title('direction of Rankine vortex field')
axis equal
axis([-xmax/2, xmax/2,-ymax/2, ymax/2])
grid on
% 
int_rv=sqrt(u.^2+v.^2);
figure
contourf(x,y,int_rv)
xlabel('x [px]')
ylabel('y [px]')
title('intensity of Rankine vortex field')
colorbar
axis equal
grid on

%% INTEGRATION OVER TIME (multiplication by time interval)
[o,r]=cart2pol(r0(:,1),r0(:,2));
uoin = (r <= R1);
uout = (r > R1);
uo = uoin+uout;
uo(uoin) =  v0*r(uoin)/R1;
uo(uout) =  v0*R1./r(uout);
uo(isnan(uo))=0;
uu = -uo.*sin(o);%by pivlab -
vv = uo.*cos(o);%by pivlab +

r_uv12=[r0 r0]+h*[-uu, -vv, uu, vv]; % particle centres [r1u r1v r2u r2v]

tic
for np=1:n_particles

   I1=I1+I0(np,nip)*exp((-(XVF-r_uv12(np,1)).^2-(YVF-r_uv12(np,2)).^2)/(1/8*dp(np,nip)^2));
   I2=I2+I0(np,nip)*exp((-(XVF-r_uv12(np,3)).^2-(YVF-r_uv12(np,4)).^2)/(1/8*dp(np,nip)^2));

end
toc

figure % check images without noise
imshow(I1, 'initialmagnification',1000)
figure
imshow(I2,'initialmagnification',1000)

I1=imcrop(I1,[extrapx+1,extrapx+1,resx-1,resy-1]); % [XMIN YMIN WIDTH HEIGHT]
I2=imcrop(I2,[extrapx+1,extrapx+1,resx-1,resy-1]); % [XMIN YMIN WIDTH HEIGHT]

I1(I1>255)=255;
I2(I2>255)=255;

I1=I1/div_f;
I2=I2/div_f;

Ifin1(nip,:,:)=I1+imnoise(I1,'gaussian',0,max(max(I1))/100);
Ifin2(nip,:,:)=I2+imnoise(I2,'gaussian',0,max(max(I2))/100);

figure % check image with noise
imshow(squeeze(Ifin1(nip,:,:)), 'initialmagnification',1000)
figure
imshow(squeeze(Ifin2(nip,:,:)),'initialmagnification',1000)

end

%% SAVING

save(name, 'Ifin1', 'Ifin2', 'ufin', 'vfin')


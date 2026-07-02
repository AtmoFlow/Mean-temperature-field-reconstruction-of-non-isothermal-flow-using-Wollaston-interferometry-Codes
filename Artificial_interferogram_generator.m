
% This script can be used to calculate a circular 2D temperature profile 
% by integrating a 3D cylindrical temperature field along the z-direction.
% This temperature profile is then used to modulate synthetic
% interferograms in x- and y-direction respectively.

clear all
clc;
close all;

%%

% data path for the input csv file containing the 3D temperature data
D = 'C:\Users\skuehne\OneDrive\Desktop\Ray\Ray tracking\Matlab Program\TEHD 3D';
S = dir(fullfile(D,'*.csv')); 
Data = cell(1, numel(S));

File1 = fullfile(S(2).folder, S(2).name);

data = readmatrix(File1);


%%

% interpolating data on new grid
Nx = 1201; Ny = 1201; Nz = 81;
xv = linspace(-10,10,Nx);
yv = linspace(-10,10,Ny);
zv = linspace(0,100,Nz);

[Xg,Yg,Zg] = meshgrid(xv,yv,zv);

F = scatteredInterpolant(data(:,1),data(:,2),data(:,3),data(:,4),'linear','nearest');
Tg = F(Xg,Yg,Zg);

% integrating along the beam path
Txy = trapz(zv, Tg, 3);
Txy = Txy / max(zv);  

Sides=size(Txy);
R=floor(Sides(1)/2);

% masking the important areas
[xp,yp]=meshgrid(-R:R,-R:R);
ind = xp.^2 + yp.^2 < (R)^2 & xp.^2 + yp.^2 > (R/2)^2;

Txy(~ind)=NaN;

Txy(all(isnan(Txy),2),:)=[];
Txy(:,all(isnan(Txy),1))=[];

Txy2=Txy-273.15;

figure
h=imagesc(flip(Txy2));
set(h, 'AlphaData', ~isnan(Txy2));
axis equal tight;
title('Artificial temperature profile'); 
clim([round(min(Txy2(:))), round(max(Txy2(:)))])
c = colorbar;
c.Ticks = linspace(round(min(Txy2(:))), round(max(Txy2(:))), 6); 
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
colormap(jet)
grid off
axis off

%%

f=50;                                 % focal length of the imaging lens [mm]
epsilon=0.26/360*2*pi;                % seperation angle of the Wollaston prism [rad]
d=round((f.*tan(epsilon)).*60);       % beam shear
lambda=532*10^-6;                     % wavelength laser [mm]
dn_dT=-3.8*10^-4;                     % thermo-optic coefficient n-T

% calculating refractive index difference and OPD
brech_index=(Tg-298.15).*dn_dT;
OPD = trapz(zv, brech_index, 3);

Sides=size(OPD);
R2=floor(Sides(1)/2);

% masking the important areas
[xp2,yp2]=meshgrid(-R2:R2,-R2:R2);
ind2 = xp2.^2 + yp2.^2 < (R2)^2 & xp2.^2 + yp2.^2 > (R2/2)^2;

OPD(~ind2)=NaN;

OPD(all(isnan(OPD),2),:)=[];
OPD(:,all(isnan(OPD),1))=[];

% calculating phase difference
phase=(((2*pi)/lambda)*OPD);

% calculating the phase gradients shearing the image with itself 
zc = NaN(size(phase,1),d);

px1=[phase,zc];
px2=[zc,phase];
dx=px1-px2;

py1=[phase; zc'];
py2=[zc'; phase];
dy=py1-py2;

% modulating the interferograms
I_x=imgaussfilt(cos(0.5.*dx).^2);
I_y=imgaussfilt(cos(0.5.*dy).^2);

I_x=filter2(fspecial('average',3),I_x);
I_y=filter2(fspecial('average',3),I_y);

% plotting phase gradients
figure
imagesc(imgaussfilt(dx));
axis equal tight;
colormap gray
title('Phase gradient x-direction'); 
colorbar

figure
imagesc(imgaussfilt(dy));
axis equal tight;
colormap gray
title('Phase gradient y-direction'); 
colorbar

% plotting synthetic interferograms
figure
h=imagesc(flip(I_x));
set(h, 'AlphaData', ~isnan(I_x));
axis equal tight;
colormap gray
title('Artificial interferogram x-gradient'); 
clim([0, 1])
c = colorbar;
c.Ticks = linspace(0, 1, 6);
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
axis off

figure
h=imagesc(flip(I_y));
set(h, 'AlphaData', ~isnan(I_y));
axis equal tight;
colormap gray
title('Artificial interferogram y-gradient'); 
clim([0, 1])
c = colorbar;
c.Ticks = linspace(0, 1, 6);
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
axis off


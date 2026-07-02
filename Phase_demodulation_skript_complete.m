
% This script can be utilized to crop two annulus pictures for the 
% X- and Y-gradient respectively into a ring shape,
% normalise them, demodulate the phase using a WFT
% algorithm, unwrap the phase, integrate them and calculate the temperature distribution.

clear all;
close all;
clc

imagefiles = dir('*.jpg');

% textbox to enter the temperature gradient and the reference temperature
prompt = {'Enter established temperature gradient:','Enter reference temperature:'};
dlgtitle = 'Input';
fieldsize = [1 45; 1 45];
definput = {'2','25'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

Gradient = str2double(answer(1));
T_ref = str2double(answer(2));

%%

% input image for the X-gradient
Pic_X = imread("2K_Hor_X2.jpg");

% input image for the Y-gradient
Pic_Y = imread("2K_Hor_Y2.jpg");

%%

% uncomment this part and plot the picture to find
% the three points for triangulation
% figure
% imshow(flip(Pic_X,1)); 

Pic_gray_X = im2double(flip(Pic_X(:,:,1),1));

% triangulation coordinates

% coordinates for example interferogram
P1_X = [1120 1065]; 
P2_X = [612 531]; 
P3_X = [1448 363];

%%

% uncomment this part and plot the picture to find
% the three points for triangulation
% figure
% imshow(flip(Pic_Y,1));

Pic_gray_Y = im2double(flip(Pic_Y(:,:,1),1));

% triangulation coordinates

% coordinates for example interferogram
P1_Y = [1035 1031]; 
P2_Y = [571 561]; 
P3_Y = [1418 341];

%%
% ring cropping for X-gradient

% triangulation for center coordiante C and radius R
C_X = NaN*[0 0]; 
C_X(1) = round( ((P2_X(1)^2 - P1_X(1)^2 + P2_X(2)^2 - P1_X(2)^2) * (P3_X(2) - P2_X(2)) ...
    - (P3_X(1)^2 - P2_X(1)^2 + P3_X(2)^2 - P2_X(2)^2) * (P2_X(2) - P1_X(2))) ...
    / (2*((P3_X(2) - P2_X(2)) * (P2_X(1) - P1_X(1)) - (P2_X(2) - P1_X(2)) * (P3_X(1) - P2_X(1)))));
C_X(2) = round( (P3_X(1)^2 - P2_X(1)^2 + P3_X(2)^2 - P2_X(2)^2 - 2*C_X(1) * (P3_X(1) - P2_X(1))) ...
    / (2 * (P3_X(2) - P2_X(2))));
 
R_X = round(((P1_X(1) - C_X(1))^2 + (P1_X(2) - C_X(2))^2)^(1/2)); % outer radius
R_Inner_X = R_X/2; % inner radius
R_Diff_X = round(R_X-R_Inner_X); % radius difference

X_north_y = C_X(1,2);      
X_north_x = C_X(1,1);

[X_nx,X_ny] = size(Pic_gray_X);
[Xind_X,Yind_X] = meshgrid(1:X_ny,1:X_nx);
x0_X=X_north_x;
y0_X=X_north_y;

for i=1:size(Yind_X,1)
    for j=1:size(Yind_X,2)
        % defining the two circular masks
        Sub1_X(i,j)=((Xind_X(i,j)-x0_X)^2+(Yind_X(i,j)-y0_X)^2)^0.5>R_X;
        Sub2_X(i,j)=((Xind_X(i,j)-x0_X)^2+(Yind_X(i,j)-y0_X)^2)^0.5<(R_X/2);

        % setting everything outside of the masks to NaN
        if Sub1_X(i,j)==1
            Pic_gray_X(i,j,:)=NaN;
        else
            Pic_gray_X(i,j,:)=Pic_gray_X(i,j,:);
        end

        if Sub2_X(i,j,:)==1
            Pic_gray_X(i,j,:)=NaN;
        else
            Pic_gray_X(i,j,:)=Pic_gray_X(i,j,:);
        end
    end
end

% deleting all rows/columns only containing NaN
Pic_gray_X(all(isnan(Pic_gray_X),2),:)=[];
Pic_gray_X(:,all(isnan(Pic_gray_X),1))=[];


% Pic_gray_X=flip(Pic_gray_X,1);
Pic_gray_X=flip(Pic_gray_X,2);

Sides_X=round(size(Pic_gray_X)/2)-1;

% define carthesian grid
[xp_X,yp_X] = meshgrid((-Sides_X(1)):(Sides_X(1)),(-Sides_X(2)):(Sides_X(2)));

figure
contourf(xp_X,yp_X,flip(Pic_gray_X,1),'LineStyle','none');   
colormap(gray);
title('X-gradient'); 
axis equal
c = colorbar;
c.Ticks = linspace(round(min(Pic_gray_X(:)),2), round(max(Pic_gray_X(:)),2), 6);
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
axis off

%%
% ring cropping for Y-gradient

% triangulation for center coordiante C and radius R
C_Y = NaN*[0 0]; 
C_Y(1) = round( ((P2_Y(1)^2 - P1_Y(1)^2 + P2_Y(2)^2 - P1_Y(2)^2) * (P3_Y(2) - P2_Y(2)) ...
    - (P3_Y(1)^2 - P2_Y(1)^2 + P3_Y(2)^2 - P2_Y(2)^2) * (P2_Y(2) - P1_Y(2))) ...
    / (2*((P3_Y(2) - P2_Y(2)) * (P2_Y(1) - P1_Y(1)) - (P2_Y(2) - P1_Y(2)) * (P3_Y(1) - P2_Y(1)))));
C_Y(2) = round( (P3_Y(1)^2 - P2_Y(1)^2 + P3_Y(2)^2 - P2_Y(2)^2 - 2*C_Y(1) * (P3_Y(1) - P2_Y(1))) ...
    / (2 * (P3_Y(2) - P2_Y(2))));
% R = NaN; 
R_Y = round(((P1_Y(1) - C_Y(1))^2 + (P1_Y(2) - C_Y(2))^2)^(1/2)); % outer radius
R_Inner_Y = R_Y/2; % inner radius
R_Diff_Y = round(R_Y-R_Inner_Y); % radius difference

Y_north_y = C_Y(1,2);      
Y_north_x = C_Y(1,1);

[nx_Y,ny_Y] = size(Pic_gray_Y);
[Xind_Y,Yind_Y] = meshgrid(1:ny_Y,1:nx_Y);
x0_Y=Y_north_x;
y0_Y=Y_north_y;

for i=1:size(Yind_Y,1)
    for j=1:size(Yind_Y,2)
        % defining the two circular masks
        Sub1_Y(i,j)=((Xind_Y(i,j)-x0_Y)^2+(Yind_Y(i,j)-y0_Y)^2)^0.5>R_Y;
        Sub2_Y(i,j)=((Xind_Y(i,j)-x0_Y)^2+(Yind_Y(i,j)-y0_Y)^2)^0.5<(R_Y/2);

        % setting everything outside of the masks to NaN
        if Sub1_Y(i,j)==1
            Pic_gray_Y(i,j,:)=NaN;
        else
            Pic_gray_Y(i,j,:)=Pic_gray_Y(i,j,:);
        end

        if Sub2_Y(i,j,:)==1
            Pic_gray_Y(i,j,:)=NaN;
        else
            Pic_gray_Y(i,j,:)=Pic_gray_Y(i,j,:);
        end
    end
end

% deleting all rows/columns only containing NaN
Pic_gray_Y(all(isnan(Pic_gray_Y),2),:)=[];
Pic_gray_Y(:,all(isnan(Pic_gray_Y),1))=[];

Sides_Y=round(size(Pic_gray_Y)/2)-1;

% define carthesian grid
[xp_Y,yp_Y] = meshgrid((-Sides_Y(1)):(Sides_Y(1)),(-Sides_Y(2)):(Sides_Y(2)));

figure
contourf(xp_Y,yp_Y,flip(flip(Pic_gray_Y,2),1),'LineStyle','none');
colormap(gray);
title('Y-gradient'); 
axis equal
c = colorbar;
c.Ticks = linspace(round(min(Pic_gray_Y(:)),2), round(max(Pic_gray_Y(:)),2), 6);
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
axis off

%%

% choose which of the two gradients should be evaluated
choice = questdlg('Which gradient will be evaluated?', 'Gradient choice', 'x', 'y','x');

switch choice
    case 'x'
        Pic_gray=Pic_gray_X;
        RR=R_X;
        
    case 'y'
        Pic_gray=Pic_gray_Y;
        RR=R_Y;
        
    otherwise
        error('No choice has been made');
end

%%

% Normalisation algorithm to increase image quality 

data=Pic_gray;
data(isnan(Pic_gray))=0;

F_data = fft2(data); F_shift_data = fftshift(F_data);
img = F_shift_data; 

[rows, cols] = size(img);
radius = 7; 

centerX = round(cols / 2);
centerY = round(rows / 2);

upperHalf = zeros(rows, cols);
rightHalf = zeros(rows, cols);

[x, y] = meshgrid(1:cols, 1:rows);
circleMask = ((x - centerX).^2 + (y - centerY).^2) <= radius^2;

upperHalf(1:centerY, :) = img(1:centerY, :);
upperHalf(circleMask) = 0;

rightHalf(:, centerX+1:end) = img(:, centerX+1:end);
rightHalf(circleMask) = 0;


F_ishift = ifftshift(upperHalf); Invx = ifft2(F_ishift);
F_anglex = angle(Invx);  
F_ishift = ifftshift(rightHalf); Invy = ifft2(F_ishift);
F_angley = angle(Invy);

irradiance=(abs(Invx).*cos(F_anglex)+abs(Invy).*cos(F_angley))./(abs(Invx)+abs(Invy));

irradiance(isnan(Pic_gray))=0;

% the normalised interferogram
figure;
colormap gray
imagesc(irradiance);
title('Normalised interferogram'); 
axis equal
axis off

%%

% WFT algorithm

f=irradiance;
f0=rescale(f);

%removing DC_offset
f=f0-mean(f0(:));

%textbox for all input parameters
prompt = {'Sigma X','Sigma Y','wxi','wxl','wxh','wyi','wyl','wyh'};
dlgtitle = 'Input WFR';
fieldsize = [1 45; 1 45; 1 45; 1 45; 1 45; 1 45; 1 45; 1 45];
definput = {'10','10','0.025','0','0.5','0.025','0','0.5'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

sigmax = str2double(answer(1));
sigmay = str2double(answer(2));
wxi = str2double(answer(3));
wxl = str2double(answer(4));
wxh = str2double(answer(5));
wyi = str2double(answer(6));
wyl = str2double(answer(7));
wyh = str2double(answer(8));

%extending kernel size
sx=round(3*sigmax);
sy=round(3*sigmay);

%defining the kernel
[y,x]=meshgrid(-sy:sy,-sx:sx);
w=1/sqrt(pi*sigmax*sigmay)*exp(-x.*x/2/sigmax/sigmax-y.*y/2/sigmay/sigmay);

%pre-definition of all arrays for phase, spatial frequencies and ridges
[m,n]=size(f); 

g_wx=zeros(m,n); 
g_wy=zeros(m,n); 
g_phase=zeros(m,n); 
g_r=zeros(m,n);

%the WFR loop
for wyt=wyl:wyi:wyh
    for wxt=wxl:wxi:wxh
        wave=w.*exp(1j*wxt*x+1j*wyt*y);

        sf=conv2(f,wave,'same');

        t=(abs(sf)>g_r);
        g_r=g_r.*(1-t)+abs(sf).*t;
        g_wx=g_wx.*(1-t)+wxt*t;
        g_wy=g_wy.*(1-t)+wyt*t;
        g_phase=g_phase.*(1-t)+angle(sf).*t;
    end
end

% re-scaling the phase
valmax=max(max(g_phase));
valmin=min(min(g_phase));
g_phase=(g_phase-valmin)/(valmax-valmin)*2*pi-pi;

g_phase(isnan(Pic_gray))=0;

% plotting the phase
figure
imagesc(g_phase)
title('Wrapped phase'); 
c = colorbar;
c.Ticks = linspace(round(min(g_phase(:))), round(max(g_phase(:))), 6); 
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
colormap(gray)
axis image
axis off

%%

% phase unwrapping

g_phase_unwrap=weighted_phase_unwrap(g_phase);

g_phase_unwrap(isnan(Pic_gray))=0;

figure
imagesc(g_phase_unwrap);
title('Unwrapped phase'); 
c = colorbar;
c.Ticks = linspace(round(min(g_phase_unwrap(:))), round(max(g_phase_unwrap(:))), 6); 
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
colormap(gray)
axis image
axis off


%%

% integration of the phase gradient

switch choice
    case 'x'
        phase_unwrap_x=g_phase_unwrap;

        [m, n] = size(phase_unwrap_x);
        phase_unwrap_X = zeros(m, n);
 
        for i = 1:m
           for j =n-1:-1:1
               phase_unwrap_X(i,j) = phase_unwrap_X(i,j+1) + phase_unwrap_x(i,j);
           end
        end

        phase_unwrap_X(isnan(Pic_gray))=NaN;

        % unwrapped phase, X-gradient
        figure; 
        imagesc(phase_unwrap_X);
        title('Integrated phase, X-Gradient');
        colormap gray
        colorbar
        axis equal
        axis off
        
    case 'y'
        phase_unwrap_y=g_phase_unwrap;

        [m, n] = size(phase_unwrap_y);
        phase_unwrap_Y = zeros(m, n);

        for i = 1:n
           for j =2:1:m
               phase_unwrap_Y(j,i) = phase_unwrap_Y(j-1,i) + phase_unwrap_y(j,i);
           end
        end
        phase_unwrap_Y(isnan(Pic_gray))=NaN;

        % unwrapped phase, Y-gradient
        figure; 
        imagesc(phase_unwrap_Y);
        title('Integrated phase, Y-Gradient');
        colormap gray
        axis equal
        axis off
        
end

%% 

switch choice
    case 'x'
        Picture_Full = phase_unwrap_X;
        
    case 'y'
        Picture_Full = phase_unwrap_Y;
        
end

Picture_Full = Picture_Full+abs(min(Picture_Full,[],'all'));

% re-scaling the phase to a temperature distribution
Picture_Full_2=Picture_Full/abs(max(Picture_Full,[],'all'))*Gradient+T_ref;

%%

% cutting out the fluid pipe
% the size of the rectangular mask needs to be chosen empirically

x = [1 RR RR 1];  % columns
y = [287 287 608 608];  % lines

mask_cut = poly2mask(x, y, size(Picture_Full_2, 1), size(Picture_Full_2, 2));

Picture_Full_2(mask_cut) = NaN;

%%

% final temperature profile
figure
h=imagesc(Picture_Full_2);
set(h, 'AlphaData', ~isnan(Picture_Full_2));
axis equal tight;
colormap jet
title('Temperature profile'); 
c = colorbar;
c.Ticks = linspace(round(min(Picture_Full_2(:)),2), round(max(Picture_Full_2(:)),2), 6);
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
axis off

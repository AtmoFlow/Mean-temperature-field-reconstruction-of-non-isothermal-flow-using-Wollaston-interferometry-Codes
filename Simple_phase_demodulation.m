
% This script can be used to generate a simple wavy interferogram 
% and subsequently demodulate it using WFT and perform a phase unwrap.

clear all;
close all;
clc;

% image size
Nx=512;
Ny=512;

[x,y]=meshgrid(0:Nx-1,0:Ny-1);

% carrier frequencies
fx=0.05;
fy=0.03;
theta=atan2(fy,fx);

A=2; % Amplitude
wx=2*pi/300; % Wave number x-direction
wy=2*pi/200; % Wave number y-direction
phi_welle=A*sin(wx*x-wy*y);

phi=2*pi*(fx*x+fy*y)+phi_welle;

I=0.5*(1+cos(phi));

figure;
imagesc(I);
title('Artificial interferogram'); 
c = colorbar;
c.Ticks = linspace(round(min(I(:))), round(max(I(:))), 6); 
c.FontSize = 20;
c.TickLabelInterpreter = 'latex';
colormap(gray)
axis image
axis off

%%

f=I;
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

%re-scaling the phase
valmax=max(max(g_phase));
valmin=min(min(g_phase));
g_phase=(g_phase-valmin)/(valmax-valmin)*2*pi-pi;

%plotting the phase
figure
imagesc(g_phase)
% title('Wrapped phase'); 
c = colorbar;
c.Ticks = linspace(round(min(g_phase(:))), round(max(g_phase(:))), 6); 
c.FontSize = 40;
c.TickLabelInterpreter = 'latex';
colormap(gray)
axis image
axis off

%%

% phase unwrapping
g_phase_unwrap=weighted_phase_unwrap(g_phase);

figure
imagesc(g_phase_unwrap);
% title('Unwrapped phase'); 
c = colorbar;
c.Ticks = linspace(round(min(g_phase_unwrap(:))), round(max(g_phase_unwrap(:))), 6); 
c.FontSize = 40;
c.TickLabelInterpreter = 'latex';
colormap(gray)
axis image
axis off


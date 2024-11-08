clear
clc
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 77e9;
N = 1024; % number of symbols
N_FFT_kx = 1024;
N_FFT_ky = 1024;
mu = 63.343e12;
fs = 9121e3;        % Sampling rate (sps)
Ts = 1/fs;          % Sampling period
z = .28; % z_target
k = 2*pi*f_0/c;
x = linspace(-0.25,0.25,N_FFT_kx) ;
y = (linspace(-1,1,N_FFT_ky))';
h = exp(-1i*2*k*sqrt(x.^2 + y.^2 + z^2));

%%
[xPointM,yPointM,~] = size(h); 
H =  fft2(h); 
%%
load('rawData3D_simple2D');  % [12, 53, 403, 256] [TX*RX Vstep Hstep N]
rawDataCal = rawData3D_simple2D(:, 1:86, 1:300 );
rawDataFFT = fft(rawDataCal,1024);
%% Range focusing to z0
tI = 4.5225e-10; % Instrument delay for range calibration (corresponds to a 6.78cm range offset)
k = round(mu*Ts*(2*z/c+tI)*1024); % corresponing range bin
sarData = squeeze(rawDataFFT(k+1,:,:));
[yPointM,xPointM] = size(sarData);
[yPointF,xPointF] = size(H);

    sarData = padarray(sarData,[0 floor((xPointF-xPointM)/2)],0,'pre');
    sarData = padarray(sarData,[0 ceil((xPointF-xPointM)/2)],0,'post');
  
    sarData = padarray(sarData,[floor((yPointF-yPointM)/2) 0],0,'pre');
    sarData = padarray(sarData,[ceil((yPointF-yPointM)/2) 0],0,'post');

sarDataFFT = fft2(sarData);
sarImage = fftshift(ifft2(sarDataFFT .* H));

xRangeT_mm = 1e-3 * (-(N_FFT_kx-1)/2 : (N_FFT_kx-1)/2); % xStepM is in mm
yRangeT_mm = 1e-3 * (-(N_FFT_ky-1)/2 : (N_FFT_ky-1)/2); % xStepM is in mm
mesh(xRangeT_mm,yRangeT_mm,abs(squeeze(sarImage)),'FaceColor','interp','LineStyle','none')

figure
imagesc(abs(squeeze(sarImage)));
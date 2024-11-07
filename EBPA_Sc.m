clear
clc
% tx_x = [0 -0.002 0]; % 3-TX
% tx_y = [0.0107 0.0147 0.0183];
tx_x = [0 0]; % 2-TX 
tx_y = [0.0117 0.0195];  % [0.0107 0.0183];
tx_y = permute(tx_y,[4,3,2,1]); 

rx_x = [0 0 0 0];
rx_y = [0.     0.0019 0.0039 0.0058];  % [0 0.0019 0.0038 0.0057];
rx_y = permute(rx_y,[4,3,1,2]);

delta_T = 0;

c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 77e9;
N = 256; % number of symbols
N0 = 2048;
N_FFT_kx = 512;
N_FFT_ky = 1024;

mu = 70.295e12;
fs = 5e6;        % Sampling rate (sps)
Ts = 1/fs;          % Sampling period
z = .25; % z_target
km = mu / c;
k = 2*pi*f_0/c;

dx = 0.98e-3;
dy = 0.98e-3;
x = dx * ([0:N_FFT_kx - 1] - N_FFT_kx / 2)'; 
y = dy * ([0:N_FFT_ky - 1] - N_FFT_ky / 2);  % 1 1 1024 1
lambda = 3e8/f_0;
rail_step_x = 0.98e-3;
rail_step_y = 7.59e-3; %8*lambda/4; 7.59e-3;
rail_step_number_x = 403;
rail_step_number_y = 53;


for ii = 0:rail_step_number_y-1
    h(:,:,ii+1,:,:) = exp(-1i*k*sqrt(x.^2 + (y-(tx_y + ii*rail_step_y)).^2 + z^2)) .* ... % tx_y : 1 3 1 1 ------>   1  3 1024 1 ---> 1 3 1024 512
        exp(-1i*k*sqrt(x.^2 + (y-(rx_y + ii*rail_step_y)).^2 + z^2)) ;     % rx_y : 4 1 1 1 ------>   4  1 1024 1 ---> 4 1 1024 512
end
H = fft(h, [], 1);
clear h

%%
load('RawDataCal.mat');  % [12, 53, 403, 256] [TX*RX Vstep Hstep N]
rawDataFFT = fft(rawDataCal,N0,4);
clear rawDataCal
%% Range focusing to z0

freq_scale = ((0:N0-1) * fs) / N0 ;% # take all positive frequencies: no fftshift
range_scale = freq_scale / (2 * km);
[o,kk] = min(abs(range_scale - z));
kk
sarData = squeeze(rawDataFFT(:,:,:,kk));
sarData = reshape (sarData,3,4,rail_step_number_y,rail_step_number_x);
sarData = sarData ([1,3],:,:,:);
sarData = permute(sarData,[4,3,1,2]); %  403    53     2     4
[yPointM,xPointM,a,b] = size(sarData);
[yPointF,xPointF,a,b,bc] = size(H);

    sarData_p = padarray(sarData,[floor((yPointF-yPointM)/2) 0],0,'pre');
    sarData_p = padarray(sarData_p,[ceil((yPointF-yPointM)/2) 0],0,'post');

S = fft(sarData_p, [], 1);
S = permute (S, [1,5,2,3,4]);
P = S.*H;
P1 = sum(sum(sum(P,3),4),5);
sarImage =fftshift( fftshift(ifft(P1, [], 1)),2);
figure
imagesc(abs(squeeze(sarImage)'));
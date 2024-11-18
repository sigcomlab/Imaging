%% EBPA method for the data wich provided by yanik in small TI MIMo radar for scissors
%% basd on the Muhammet Emin Yanik paper 'Development and Demonstration of MIMO-SAR mmWave Imaging Testbeds'.
clear
clc
load('rawData3D_simple2D');  % % [512, 100, 403] [N Vstep Hstep]
rawDataCal = permute(rawData3D_simple2D, [4,2,3,1]);
rawDataCal = rawDataCal(:, :, :, 1:256);
%% TX RX position
% tx_x = [0 -0.002 0]; % 3-TX
% tx_y = [0.0107 0.0147 0.0183];
% tx_x = [0 0]; % 2-TX 
% tx_y = [0.0117 0.0195]; 
tx_x = [0]; % 2-TX 
tx_y = [0];  %[0.0117]; 
tx_y = permute(tx_y,[4,3,2,1]); 
% 
% rx_x = [0 0 0 0];
% rx_y = [0  0.0019 0.0039 0.0058]; 
rx_x = [0];
rx_y = [0]; 
rx_y = permute(rx_y,[4,3,1,2]);

delta_T = 0; % The distance between TXs and RXs
%% radar properties
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 77e9; % start frequency
N = size(rawDataCal, 4); % number of symbols
N0 = 2048; 
N_FFT_kx = 512; % number of symbols in x-axis
N_FFT_ky = 1024; % number of symbols in y-axis

mu = 63.343e12;  % 70.295e12; % slope
fs = 9121e3;   % 5e6;        % Sampling rate (sps)
Ts = 1/fs;          % Sampling period
z = .25; % z_target distance
km = mu / c;
k = 2*pi*f_0/c;
%% Rail properties
dx = 0.98e-3;  % each step in x-axis on the rail
dy = 0.98e-3;  % each step in y-axis on the rail
x = dx * ([0:N_FFT_kx - 1] - N_FFT_kx / 2)';  % steps in x-axis 
y = dy * ([0:N_FFT_ky - 1] - N_FFT_ky / 2);  % 1 1 1024 1
lambda = 3e8/f_0;
rail_step_x = 0.98e-3;
rail_step_y = rail_step_x;  % 7.59e-3; % 8*lambda/4; 
%% loading the input data
% load('RawDataCal.mat');  % [12, 53, 403, 256] [TX*RX Vstep Hstep N]

% rawDataCal = rawDataCal(:,25:36,:,:); % by reducing the number of vertical % steps, it does not work
rail_step_number_x = size(rawDataCal,3);
rail_step_number_y = size(rawDataCal,2);
%% Making H as match filter equation
for ii = 0:rail_step_number_y-1
    h(:,:,ii+1,:,:) = exp(-1i*k*sqrt(x.^2 + (y-(tx_y + ii*rail_step_y)).^2 + z^2)) .* ... 
        exp(-1i*k*sqrt(x.^2 + (y-(rx_y + ii*rail_step_y)).^2 + z^2)) ;    
end
H = fft(h, [], 1); % [512 1024 53 2 4]
clear h

%% making S as FFT form of imported signal
rawDataFFT = fft(rawDataCal,N0,4); % [12, 53, 403, 2048] [TX*RX Vstep Hstep N0]
clear rawDataCal
%% Range focusing to z0
freq_scale = ((0:N0-1) * fs) / N0 ;% # take all positive frequencies: no fftshift
range_scale = freq_scale / (2 * km);
[o,kk] = min(abs(range_scale - z));
%% Arrangement of imported data
sarData = squeeze(rawDataFFT(:,:,:,kk)); % [12, 53, 403, 1] [TX*RX Vstep Hstep kk] selecting one N corresponding the FFT_Range
sarData = reshape (sarData,1,1,rail_step_number_y,rail_step_number_x);  % [3, 4, 53, 403] 
sarData = sarData ([1],:,:,:); % [2, 4, 53, 403] 
sarData = permute(sarData,[4,3,1,2]); % [403 53 2 4]
[yPointM,xPointM,a,b] = size(sarData);
[yPointF,xPointF,a,b,bc] = size(H);
%% Apply padding
sarData_p = padarray(sarData,[floor((yPointF-yPointM)/2) 0],0,'pre');
sarData_p = padarray(sarData_p,[ceil((yPointF-yPointM)/2) 0],0,'post'); % [512 53 2 4]

S = fft(sarData_p, [], 1); % [512 53 2 4]
S = permute (S, [1,5,2,3,4]); % [512 1 53 2 4]
P = S.*H; % [512 1024 53 2 4]
P1 = sum(sum(sum(P,3),4),5); % [512 1024]
sarImage =fftshift( fftshift(ifft(P1, [], 1)),2); % [512 1024]
figure
imagesc(abs(squeeze(sarImage)'));

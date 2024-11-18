%% EBPA method for the data wich provided in big TI MIMo radar for Pliers/scissors based on the method in article 'Development and Demonstration of MIMO-SAR mmWave Imaging Testbeds'
clear
clc
%% Selection the data between pliers or scissors
% load('scissors.mat') % [300 12 16 512], [x-step-on-rail TX RX N]
% z = 0.4 ; % z_target for scissors

load('pliers_calibrated.mat');  % [300 12 16 512] [Vstep_on_rail TX RX N]
z = 0.75; % z_target for pliers (The distance between target and radar)
%% TX RX position
tx_x = [0	0	0	0	0	0	0	0	0	0.00191082802547771	0.00764331210191083	0.0114649681528662];
tx_y = [0.0152866242038217	0.0229299363057325	0.0305732484076433	0.0382165605095541	0.0458598726114650	0.0535031847133758	0.0611464968152866	0.0687898089171975	0.0764331210191083	0.0324840764331210	0.0343949044585987	0.0363057324840764];
rx_x = [0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420];
rx_y = [0	0.00191082802547771	0.00382165605095541	0.00573248407643312	0.0210191082802548	0.0229299363057325	0.0248407643312102	0.0267515923566879	0.0878980891719745	0.0898089171974522	0.0917197452229299	0.0936305732484077	0.0955414012738854	0.0974522292993631	0.0993630573248408	0.101273885350318];
% load('tx_x.mat'),load("tx_y.mat"),load("rx_x.mat"),load("rx_y.mat"),
% Removing non-aligned TXs
tx_x (10:12) = []; 
tx_y (10:12) = [];
% Appropriate allocation of TX and RX on the 4-D matrix
tx_y = permute(tx_y,[4,3,2,1]) + randn(1,1,9)*0.000001;  % [1 1 9] --> [1 1 TX]
rx_y = permute(rx_y,[4,3,1,2]) + randn(1,1,1,16)*0.000001; % [1 1 1 16] --> [1 1 1 RX]
delta_T = 0; % Distance between TX and RX
%% Radar properties
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 78.5e9;
N = 512; % number of symbols
N0 = 2048;
N_FFT_kx = 512;  % number of symbols in x-axis
N_FFT_ky = 1024; % number of symbols in y-axis
mu = 86e12; % Slope
fs = 10e6;        % Sampling rate (sps)
Ts = 1/fs;          % Sampling period
km = mu / c;
k = 2*pi*f_0/c;
% rail steps
dx = 0.98e-3; % each step in x-axis on the rail
dy = 0.975e-3; % each step in y-axis on the rail
% target domain definition
x = dx * ((0:N_FFT_kx - 1) - N_FFT_kx / 2)';   % [512 1] --> [N_FFT_kx 1]
y = dy * ((0:N_FFT_ky - 1) - N_FFT_ky / 2);  % [1 1024] --> [1 N_FFT_ky]
lambda = 3e8/f_0;
rail_step_x = 0.975e-3;
rail_step_y = 7.59e-3; %8*lambda/4; 7.59e-3;
rail_step_number_x = 300;
rail_step_number_y = 1;
%% Matching Filter definition
for ii = 0:rail_step_number_y-1
    h(:,:,ii+1,:,:) = exp(-1i*k*sqrt(x.^2 + (y-(tx_y + ii*rail_step_y)).^2 + z^2)) .* ... % tx_y : 1 3 1 1 ------>   1  3 1024 1 ---> 1 3 1024 512
        exp(-1i*k*sqrt(x.^2 + (y-(rx_y + ii*rail_step_y)).^2 + z^2)) ;     % rx_y : 4 1 1 1 ------>   4  1 1024 1 ---> 4 1 1024 512
end
H = fft(h, [], 1); % [512 1024 1 9 16] 
clear h
%% FFT of the data (pliers, scissors)
rawDataFFT = fft(s_full,N0,4); % [300 12 16 2048] [Vstep_on_rail TX RX N0]
clear rawDataCal
%% Range focusing to z0
freq_scale = ((0:N0-1) * fs) / N0 ;% # take all positive frequencies: no fftshift
range_scale = freq_scale / (2 * km);
[o,kk] = min(abs(range_scale - z));
% kk
%% Data arrangment
sarData = squeeze(rawDataFFT(:,:,:,kk)); % [300 12 16 1] [Vstep_on_rail TX RX kk] after selection of max rang in range FFT
sarData = sarData (:,[1:9],:); % [300 9 16] [Vstep_on_rail TX RX]
sarData = permute(sarData,[1,4,5,2,3]); % [300 1 1 9 16] [Vstep_on_rail TX RX N]
[yPointM,xPointM,a,b] = size(sarData);
[yPointF,xPointF,a,b,bc] = size(H);
%% Array padding
sarData_p = padarray(sarData,[floor((yPointF-yPointM)/2) 0],0,'pre');
sarData_p = padarray(sarData_p,[ceil((yPointF-yPointM)/2) 0],0,'post'); % [512 1 1 9 16] [Vstep_on_rail 1 1 TX RX]
%% Applying FFT on main signal
S = fft(sarData_p, [], 1); % [512 1 1 9 16] [Vstep_on_rail TX RX N]
P = S.*H; % [512 1024 1 9 16] 
P1 = sum(sum(sum(P,3),4),5); % [512 1024]
sarImage =fftshift( fftshift(ifft(P1, [], 1)),2); % [512 1024]
% figure
sarImage = flip(sarImage,2);  % [512 1024]
imagesc(abs(squeeze(sarImage)'));

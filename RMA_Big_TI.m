%% RMA method based on the data of Big TI in our Lab  based on the method in 'Development and Demonstration of MIMO-SAR mmWave Imaging Testbeds'.
clear
clc
close all
load('pliers_calibrated.mat') % [300 12 16 512], [x-step-on-rail TX RX N]
% load('tx_x.mat'),load("tx_y.mat"),load("rx_x.mat"),load("rx_y.mat"),
tx_x = [0	0	0	0	0	0	0	0	0	0.00191082802547771	0.00764331210191083	0.0114649681528662];
tx_y = [0.0152866242038217	0.0229299363057325	0.0305732484076433	0.0382165605095541	0.0458598726114650	0.0535031847133758	0.0611464968152866	0.0687898089171975	0.0764331210191083	0.0324840764331210	0.0343949044585987	0.0363057324840764];
rx_x = [0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420	0.0649681528662420];
rx_y = [0	0.00191082802547771	0.00382165605095541	0.00573248407643312	0.0210191082802548	0.0229299363057325	0.0248407643312102	0.0267515923566879	0.0878980891719745	0.0898089171974522	0.0917197452229299	0.0936305732484077	0.0955414012738854	0.0974522292993631	0.0993630573248408	0.101273885350318];
tx_x (10:12) = []; 
tx_y (10:12) = [];
s = s_full;
s = s(:,[1:9],:,:);  % Removing not align TX ones
delta_T = rx_x (1);
plot(tx_x, tx_y  ,'O'), hold on, grid on
plot(rx_x, rx_y ,'*')

%% Making VAs and narrow down the non-overlapped ones
i = 0;
for x = 1:length (tx_x)
    for y = 1: length(rx_x)
        VA_x(i+1) = (tx_x(x) + rx_x(y)) / 2;
        VA_y(i+1) = (tx_y(x) + rx_y(y)) / 2;
        plot(VA_x(i+1), VA_y(i+1) ,'.')
        i = i + 1;
        C(:,i) = [x y];
    end
end
[UniqX, iX] = unique(round(VA_y,4));
VA_useful =  C(:,iX); % [TX_index  RX_index] 
TX_useful = VA_useful(1,:); % index of usefule Tx for making 86 VAs [1*86]
Rx_usful = VA_useful(2,:);  % index of usefule Rx for making 86 VAs [1*86]
%% Define Frequency Spectrum
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 78.5e9;
N = 512; % number of symbols
N_FFT_kx = 512;  % number of symbols in x-axis
N_FFT_ky = 1024; % number of symbols in y-axis
K = 86e12;
f_s = 10e6; % sampeling frequency
f = f_0 + (0:N-1)*K/f_s;  
K = 2*pi*f/c;
k = reshape(K,1,1,[]);

%% kx,kx,kz
kx = linspace(-(pi/1e-3),(pi/1e-3),N_FFT_kx); % kX-Domain
ky = (linspace(-(pi/1e-3),(pi/1e-3),N_FFT_ky)).'; % kY-Domain

kz = single(sqrt((2*k).^2 - kx.^2 - ky.^2));
kz((kx.^2 + ky.^2) > (2*k).^2) = 0;

%% calculation the distances
x = 0.15;  % x_target;
y = 0.0; % y_target; 
z = .835; % z_target;

rail_step_x = 0.975e-3; % each step in x-axis on the rail
rail_step_y = 0.975e-3*86;  % each step in y-axis on the rail
rail_step_number_x = 300; % number of steps in x-axis
rail_step_number_y = 1; % number of steps in y-axis

rail_x = rail_step_x * (1:rail_step_number_x); % steps in x-axis
rail_y = rail_step_y * (1:rail_step_number_y); % steps in y-axis
% plot(rail_x,zeros(1,length(rail_x)), '.', 'linewidth',2)
% plot(zeros(1,length(rail_y)),rail_y, '.', 'linewidth',2)
RawDataCal = s;

for ii = 0:rail_step_number_y-1
    for i = 0:rail_step_number_x-1
        for j = 1:length(tx_x)
            for l = 1:length(rx_x)
                    x_prime = (tx_x(j) + rx_x(l)) / 2 + i * rail_step_x;
                    y_prime = (tx_y(j) + rx_y(l)) / 2 + ii* rail_step_y;
                    % plot(x_prime,y_prime,'o','LineWidth',2), hold on, grid on

                    R_T = sqrt ((x-(x_prime + delta_T/2))^2 + (y-tx_y(j))^2 + (z)^2);
                    R_R = sqrt ((x-(x_prime - delta_T/2))^2 + (y-rx_y(l))^2 + (z)^2);
                    R = sqrt ((x-x_prime)^2 + (y-y_prime)^2 + (z)^2);
                    for k = 1:length(K)
                        s_hat_0_mono = exp(-1i*2*K(k)*R);
                        s_hat_0_multi = exp(-1i*K(k)* (R_T+R_R));
                        s_tilde(i+1,j,l,k) = RawDataCal(i+1,j,l,k) * (s_hat_0_mono/s_hat_0_multi);
                    end
            end
        end
    end
end
% [300 12 16 512], [x-step-on-rail TX RX N]
% tx_idx = [0 0 0 0 1 1 1 1 2 2 2 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 ...
%     8 8 8 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8];
% rx_idx = [ 0  1  2  3  0  1  2  3  0  1  2  4  5  6  7  4  5  6  7  4  5  6  7  4 ...
%     5  6  7  4  5  6  7  4  5  6  7  4  5  6  7  4  5  6  7  4  5  6  8  9 10 11 12 13 ...
%     14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15];
% for idx = 1:86
%     s_tilde_uniform(:, idx, :) = s_tilde(:,tx_idx(idx)+1, rx_idx(idx)+1, :);   % [300 86*1 512], [x-step on rail TX RX N] after removing overlapped ones
% end
for idx = 1:length(TX_useful)
    s_tilde_uniform(:, idx, :) = s_tilde(:,TX_useful(idx), Rx_usful(idx), :);   % [300 86*1 512], [x-step on rail TX RX N] after removing overlapped ones
end
% s_tilde_uniform = reshape (s_tilde,9*16,rail_step_number_x,N);

[xPointM,yPointM,~] = size(s_tilde_uniform);

s_tilde_uniform_padd = padarray(s_tilde_uniform,[floor((N_FFT_kx-xPointM)/2) 0],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[ceil((N_FFT_kx-xPointM)/2) 0],0,'post');

s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 floor((N_FFT_ky-yPointM)/2)],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 ceil((N_FFT_ky-yPointM)/2)],0,'post'); % [512 1024 512]

S_tilde_kx_ky = fftshift(fftshift(fft2(s_tilde_uniform_padd),1),2); % [512 1024 512]

phaseFactor = exp(-1i*z*kz); 
phaseFactor= permute(phaseFactor,[2 1 3]);  % [512 1024 512]
S_tilde_kx_ky = S_tilde_kx_ky .* phaseFactor;  % [512 1024 512]
S_tilde_kx_ky = sum(S_tilde_kx_ky,3);  % [512 1024]
sarImage = ifft2(S_tilde_kx_ky);  % [512 1024]
sarImage = flip(sarImage,2);  % [512 1024]

xRangeT_mm = 1e-3 * (-(N_FFT_kx-1)/2 : (N_FFT_kx-1)/2); % xStepM is in mm
yRangeT_mm = 1e-3 * (-(N_FFT_ky-1)/2 : (N_FFT_ky-1)/2); % xStepM is in mm

% xRangeT_mm = xRangeT_mm(indXpartT);
% yRangeT_mm = yRangeT_mm(indYpartT);
% mesh(xRangeT_mm,yRangeT_mm,abs(squeeze(sarImage)),'FaceColor','interp','LineStyle','none')
imagesc(abs(squeeze(sarImage)))
view(2)

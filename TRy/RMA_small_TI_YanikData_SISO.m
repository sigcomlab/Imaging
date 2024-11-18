%% RMA method for the data wich provided by yanik as shapes in small TI MIMo radar
%% basd on the Muhammet Emin Yanik paper 'Development and Demonstration of MIMO-SAR mmWave Imaging Testbeds'.
% This cod is not completed
clear
clc
%% TX RX position
tx_x = [0 -0.002 0];  % 3-TX
tx_y = [0.0107 0.0147 0.0186]; 
which_tx = 2;
which_rx = 4;

rx_x = [0 0 0 0];
rx_y = [0 0.0019 0.0038 0.0057];
%% Importing the data and rail properties
load('rawData3D_simple2D.mat'); % [512 100 403] [TX*RX vertical-step horizental-step N]
[N,rail_step_number_y,rail_step_number_x,~] = size(rawData3D_simple2D);
rawDataCal = rawData3D_simple2D;

rail_step_x = 1e-3; % each step in x-axis on the rail
rail_step_y = 1e-3; % each step in y-axis on the rail
%% Define Frequency Spectrum
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 77e9;
N_FFT_kx = 1024;
N_FFT_ky = 1024;
K = 63.343e12;
f_s = 9121e3; 
Ts = 1/f_s;
f = f_0 + (0:N-1)*K/f_s; 
K = 2*pi*f/c;
k = reshape(K,1,1,[]);

%% kx,kx,kz
kx = linspace(-(pi/1e-3),(pi/1e-3),N_FFT_kx); % kX-Domain

ky = (linspace(-(pi/1e-3),(pi/1e-3),N_FFT_ky)).'; % kY-Domain

kz = single(sqrt((2*k).^2 - kx.^2 - ky.^2));
kz((kx.^2 + ky.^2) > (2*k).^2) = 0;

%% RMA Method for imaging
x = 0.;  % 19747;  % x_target
y = 0.;  % 2011; % y_target
z = 0.28; % z_target

RawDataCal = permute (rawDataCal,[5,4,2,3,1]); % [1 1 100 403 512] [TX RX vertical-step horizental-step N]
for ii = 0:rail_step_number_y-1
    for i = 0:rail_step_number_x-1
        for j = 1:length(tx_x(which_tx))
            for l = 1:length(rx_x(which_rx))
                    R_T = sqrt ((x-(i * rail_step_x))^2 + (y-(tx_y(j)+ii*rail_step_y))^2 + (z)^2);
                    R_R = sqrt ((x-(i * rail_step_x))^2 + (y-(rx_y(l)+ii*rail_step_y))^2 + (z)^2);
                    x_prime = (tx_x(j)+rx_x(l))/2;
                    y_prime = (tx_y(j)+rx_y(l))/2;
                    R = sqrt ((x-(x_prime + i * rail_step_x))^2 + (y-(y_prime+ii*rail_step_y))^2 + (z)^2);
                    for k = 1:length(K)
                        s_hat_0_mono = exp(-1i*2*K(k)*R);
                        s_hat_0_multi = exp(-1i*K(k)* (R_T+R_R));
                        s_tilde(j,l,ii+1,i+1,k) = RawDataCal(j,l,ii+1,i+1,k) * (s_hat_0_mono/s_hat_0_multi);
                    end
            end
        end
    end
end
%% Arrangement of imported data
s_tilde = s_tilde([1],:,:,:,:); % [2 4 53 403 256] [TX RX vertical-step horizental-step N]
s_tilde = reshape (s_tilde,1,[],rail_step_number_x,N); % [2*4 53 403 256] [TX*RX vertical-step horizental-step N]
s_tilde_uniform = reshape (s_tilde,[],rail_step_number_x,N);
[xPointM,yPointM,~] = size(s_tilde_uniform);
%% Apply padding
s_tilde_uniform_padd = padarray(s_tilde_uniform,[floor((N_FFT_kx-xPointM)/2) 0],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[ceil((N_FFT_kx-xPointM)/2) 0],0,'post');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 floor((N_FFT_ky-yPointM)/2)],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 ceil((N_FFT_ky-yPointM)/2)],0,'post'); % [512 512 256]

S_tilde_kx_ky = fftshift(fftshift(fft2(s_tilde_uniform_padd),1),2); % [512 512 256]
clear s_tilde_uniform_padd,
%% Apply phase factor
phaseFactor = exp(-1i*z*kz); % [512 512 256]
phaseFactor= permute(phaseFactor,[2 1 3]); % [512 512 256]
S_tilde_kx_ky = S_tilde_kx_ky .* phaseFactor; % [512 512 256]
S_tilde_kx_ky = sum(S_tilde_kx_ky,3); % [512 512]
sarImage = ifft2(S_tilde_kx_ky); % [512 512]
sarImage = flip(sarImage,2);

xRangeT_mm = 1e-3 * (-(N_FFT_kx-1)/2 : (N_FFT_kx-1)/2); % xStepM is in mm
yRangeT_mm = 1e-3 * (-(N_FFT_ky-1)/2 : (N_FFT_ky-1)/2); % xStepM is in mm

% xRangeT_mm = xRangeT_mm(indXpartT);
% yRangeT_mm = yRangeT_mm(indYpartT);
 mesh(xRangeT_mm,yRangeT_mm,abs(squeeze(sarImage)),'FaceColor','interp','LineStyle','none')
% imagesc(abs(squeeze(sarImage)))
view(2)

clear
clc
% close all
load('RawDataCal.mat');
tx_x = [0 -0.002 0]; % 3-TX
% tx_y = [0.0107 0.0147 0.0183];
tx_y = [0.0117 0.0147 0195];

% tx_x = [0 0]; % 2-TX
% tx_y = [0.0107 0.0183];

rx_x = [0 0 0 0];
rx_y = [0 0.0019 0.0038 0.0057];
delta_T = 0;
% plot(tx_x, tx_y  ,'O'), hold on, grid on
% plot(rx_x, rx_y ,'*')

%% Making VAs
i = 0;
for x = 1:length (tx_x)
    for y = 1: length(rx_x)
        VA_x(i+1) = (tx_x(x) + rx_x(y)) / 2;
        VA_y(i+1) = (tx_y(x) + rx_y(y)) / 2;
        i = i + 1;
    end
end
% plot(VA_x, VA_y ,'+')

%% Define Frequency Spectrum
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 77e9;
N = 256; % number of symbols
N_FFT_kx = 512;
N_FFT_ky = 512;
K = 70.295e12;
f_s = 5e6;
f = f_0 + (0:N-1)*K/f_s; 
K = 2*pi*f/c;
k = reshape(K,1,1,[]);
%% kx,kx,kz
kx = linspace(-(pi/1e-3),(pi/1e-3),N_FFT_kx); % kX-Domain

ky = (linspace(-(pi/1e-3),(pi/1e-3),N_FFT_ky)).'; % kY-Domain

kz = single(sqrt((2*k).^2 - kx.^2 - ky.^2));
kz((kx.^2 + ky.^2) > (2*k).^2) = 0;

%% calculation the distances
x = 0.19747;  % x_target; % i will define it
y = 0.2011; % y_target; % i will define it
z = 0.25; % z_target; % i will define it

rail_step_x = 0.98e-3;
rail_step_y = 7.59e-3; %7.59e-3;
rail_step_number_x = 403;
rail_step_number_y = 53;

rail_x = rail_step_x * [1:rail_step_number_x];
rail_y = rail_step_y * [1:rail_step_number_y];
plot(rail_x,zeros(1,length(rail_x)), '.', 'linewidth',2)
plot(zeros(1,length(rail_y)),rail_y, '.', 'linewidth',2)


% x_prime = rx_x (1)/2;
% y_prime = 0.0059;

RawDataCal = reshape (rawDataCal,length(tx_x),length(rx_x),[],rail_step_number_x,N);
% RawDataCal = permute(RawDataCal , [2 1 3 4 5]);
for ii = 0:rail_step_number_y-1
    for i = 0:rail_step_number_x-1
        for j = 1:length(tx_x)
            for l = 1:length(rx_x)
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

% s_tilde: [3,4,53,403,256]
% s_tilde: [12,53,403,256]
% s_tilde_uniform: [8*53,403,256]
% s_tilde_uniform_padd: [512,512,256]
s_tilde = s_tilde([1,3],:,:,:,:);
s_tilde = reshape (s_tilde,8,[],403,256);
s_tilde_uniform = reshape (s_tilde,[],403,256);


[xPointM,yPointM,~] = size(s_tilde_uniform);

s_tilde_uniform_padd = padarray(s_tilde_uniform,[floor((N_FFT_kx-xPointM)/2) 0],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[ceil((N_FFT_kx-xPointM)/2) 0],0,'post');

s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 floor((N_FFT_ky-yPointM)/2)],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 ceil((N_FFT_ky-yPointM)/2)],0,'post');

S_tilde_kx_ky = fftshift(fftshift(fft2(s_tilde_uniform_padd),1),2);
clear s_tilde_uniform_padd,

phaseFactor = exp(-1i*z*kz);
phaseFactor= permute(phaseFactor,[2 1 3]);
S_tilde_kx_ky = S_tilde_kx_ky .* phaseFactor;
S_tilde_kx_ky = sum(S_tilde_kx_ky,3);
sarImage = ifft2(S_tilde_kx_ky);
sarImage = flip(sarImage,2);

xRangeT_mm = 1e-3 * (-(N_FFT_kx-1)/2 : (N_FFT_kx-1)/2); % xStepM is in mm
yRangeT_mm = 1e-3 * (-(N_FFT_ky-1)/2 : (N_FFT_ky-1)/2); % xStepM is in mm

% xRangeT_mm = xRangeT_mm(indXpartT);
% yRangeT_mm = yRangeT_mm(indYpartT);
 mesh(xRangeT_mm,yRangeT_mm,abs(squeeze(sarImage)),'FaceColor','interp','LineStyle','none')
% imagesc(abs(squeeze(sarImage)))
view(2)
clear
clc
% close all
% load('/media/lab/Data/matlab code Yanik_Symeo/codes_Development/pliers_calibrated.mat')
load('pliers_calibrated.mat')
load('tx_x.mat'),load("tx_y.mat"),load("rx_x.mat"),load("rx_y.mat"),
% load('s.mat');
tx_x (10:12) = []; 
tx_y (10:12) = [];
s = s_full;
% s = permute(data,[1,3,4,2]);
% s(:,:,:,:,2) = [];
s = s(:,[1:9],:,:); 
% s =  permute(s,[1 3 4 2]); %  [300    12    16   512]
delta_T = rx_x (1);
x_prime = rx_x (1)/2;
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
f_0 = 78.5e9;
N = 512; % number of symbols
N_FFT_kx = 512;
N_FFT_ky = 1024;
K = 86e12;
f_s = 10e6;
f = f_0 + (0:N-1)*K/f_s; 
K = 2*pi*f/c;
k = reshape(K,1,1,[]);

%% kx,kx,kz
kx = linspace(-(pi/1e-3),(pi/1e-3),N_FFT_kx); % kX-Domain
ky = (linspace(-(pi/1e-3),(pi/1e-3),N_FFT_ky)).'; % kY-Domain

kz = single(sqrt((2*k).^2 - kx.^2 - ky.^2));
kz((kx.^2 + ky.^2) > (2*k).^2) = 0;

%% calculation the distances
x = 0.15; % 0.19747;  % x_target; % i will define it
y = 0.0; % y_target; % i will define it
z = .835; % z_target; % i will define it

rail_step_x = 0.975e-3;
rail_step_y = 0.975e-3*86;  %7.59e-3; %7.59e-3;
rail_step_number_x = 300;
rail_step_number_y = 1;

rail_x = rail_step_x * [1:rail_step_number_x];
rail_y = rail_step_y * [1:rail_step_number_y];
% plot(rail_x,zeros(1,length(rail_x)), '.', 'linewidth',2)
% plot(zeros(1,length(rail_y)),rail_y, '.', 'linewidth',2)
RawDataCal = s;

% RawDataCal = reshape (rawDataCal,length(tx_x),length(rx_x),[],rail_step_number_x,N);
% RawDataCal = permute(RawDataCal , [2 1 3 4 5]);
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
                    % R_T = sqrt ((x-(i * rail_step_x))^2 + (y-(tx_y(j)+ii*rail_step_y))^2 + (z)^2);
                    % R_R = sqrt ((x-(i * rail_step_x))^2 + (y-(rx_y(l)+ii*rail_step_y))^2 + (z)^2);
                    % R = sqrt ((x-(x_prime + i * rail_step_x))^2 + (y-(y_prime+ii*rail_step_y))^2 + (z)^2);
                    for k = 1:length(K)
                        s_hat_0_mono = exp(-1i*2*K(k)*R);
                        s_hat_0_multi = exp(-1i*K(k)* (R_T+R_R));
                        s_tilde(i+1,j,l,k) = RawDataCal(i+1,j,l,k) * (s_hat_0_mono/s_hat_0_multi);
                    end
            end
        end
    end
end

%: [9,16,1,300,256]
%  [144,1,300,256]
%  [144*1,300,256]
% [512,512,256]

% s_tilde = permute (s_tilde,[2,3,1,4]);
tx_idx = [0 0 0 0 1 1 1 1 2 2 2 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 ...
    8 8 8 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8];
rx_idx = [ 0  1  2  3  0  1  2  3  0  1  2  4  5  6  7  4  5  6  7  4  5  6  7  4 ...
    5  6  7  4  5  6  7  4  5  6  7  4  5  6  7  4  5  6  7  4  5  6  8  9 10 11 12 13 ...
    14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15 12 13 14 15];
for idx = 1:86
    s_tilde_uniform(:, idx, :) = s_tilde(:,tx_idx(idx)+1, rx_idx(idx)+1, :);
end


% s_tilde_uniform = reshape (s_tilde,9*16,rail_step_number_x,N);


[xPointM,yPointM,~] = size(s_tilde_uniform);

s_tilde_uniform_padd = padarray(s_tilde_uniform,[floor((N_FFT_kx-xPointM)/2) 0],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[ceil((N_FFT_kx-xPointM)/2) 0],0,'post');

s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 floor((N_FFT_ky-yPointM)/2)],0,'pre');
s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 ceil((N_FFT_ky-yPointM)/2)],0,'post');

S_tilde_kx_ky = fftshift(fftshift(fft2(s_tilde_uniform_padd),1),2);
% S_tilde_kx_ky =fft2(s_tilde_uniform_padd);

% clear s_tilde_uniform_padd,

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
% mesh(xRangeT_mm,yRangeT_mm,abs(squeeze(sarImage)),'FaceColor','interp','LineStyle','none')
imagesc(abs(squeeze(sarImage)))
view(2)


%% Small TI
% c = 299792458; % physconst('lightspeed'); in m/s
% f_0 = 76e9;
% N = 512; % number of symbols
% N_FFT_kx = 512;
% N_FFT_ky = 1024;
% K = 87e12;
% f_s = 10e6;
% f = f_0 + (0:N-1)*K/f_s; 
% K = 2*pi*f/c;
% k = reshape(K,1,1,[]);
% %% kx,kx,kz
% kx = linspace(-(pi/1e3),(pi/1e-3),N_FFT_kx); % kX-Domain
% 
% ky = (linspace(-(pi/1e3),(pi/1e3),N_FFT_ky)).'; % kY-Domain
% 
% kz = single(sqrt((2*k).^2 - kx.^2 - ky.^2));
% % kz((kx.^2 + ky.^2) > (2*k).^2) = 0;
% 
% %% calculation the distances
% x = 0.1462; % x_target; % i will define it
% y = 0.0400; % y_target; % i will define it
% z = 0.9; % z_target; % i will define it
% rail_step = 0.975e-3;
% rail_step_number = 300;
% rail_x = rail_step * [1:rail_step_number];
% rail_y = zeros(1,length(rail_x));
% plot(rail_x,rail_y,'.')
% y_prime = 0.0506;
% 
% for i = 0:rail_step_number-1
%     for j = 1:length(tx_x)
%         for l = 1:length(rx_x)
%             for k = 1:length(K)
%                 R_T = sqrt ((x-((x_prime + i * rail_step) + delta_T/2))^2 + (y-tx_y(j))^2 + (z)^2);
%                 R_R = sqrt ((x-((x_prime + i * rail_step) - delta_T/2))^2 + (y-rx_y(l))^2 + (z)^2);
%                 R = sqrt ((x-(x_prime + i * rail_step))^2 + (y-y_prime)^2 + (z)^2);
%                 % s(i+1,j,l,k) = exp(-1i*k*R_T) * exp(-1i*k*R_R) / (R_T * R_R); 
%                 s_hat_0_mono = exp(-1i*2*K(k)*R);
%                 s_hat_0_multi = exp(-1i*K(k)* (R_T+R_R));
%                 s_tilde(i+1,j,l,k) = s(i+1,j,l,k) * (s_hat_0_mono/s_hat_0_multi);
%             end
%         end
%     end
% end
% % s_tilde: [300,12,16,512]
% % s_tilde_uniform: [300,9*16=144,512]
% % s_tilde_uniform_padd: [512,512,512]
% 
% s_tilde_uniform = reshape(s_tilde(:,[1:9],:,:),rail_step_number,[],N);
% 
% [xPointM,yPointM,~] = size(s_tilde_uniform);
% 
% s_tilde_uniform_padd = padarray(s_tilde_uniform,[floor((N_FFT_kx-xPointM)/2) 0],0,'pre');
% s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[ceil((N_FFT_kx-xPointM)/2) 0],0,'post');
% 
% s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 floor((N_FFT_ky-yPointM)/2)],0,'pre');
% s_tilde_uniform_padd = padarray(s_tilde_uniform_padd,[0 ceil((N_FFT_ky-yPointM)/2)],0,'post');
% 
% S_tilde_kx_ky = fftshift(fftshift(fft2(s_tilde_uniform_padd),1),2);
% clear s_tilde_uniform_padd,
% 
% phaseFactor = exp(-1i*z*kz);
% S_tilde_kx_ky = S_tilde_kx_ky .* phaseFactor;
% S_tilde_kx_ky = sum(S_tilde_kx_ky,3);
% sarImage = ifft2(S_tilde_kx_ky);
% % sarImage = flip(sarImage,2);
% 
% xRangeT_mm = 1e-3 * (-(N_FFT_kx-1)/2 : (N_FFT_kx-1)/2); % xStepM is in mm
% yRangeT_mm = 1e-3 * (-(N_FFT_ky-1)/2 : (N_FFT_ky-1)/2); % xStepM is in mm
% 
% % xRangeT_mm = xRangeT_mm(indXpartT);
% % yRangeT_mm = yRangeT_mm(indYpartT);
% mesh(xRangeT_mm,yRangeT_mm,abs(squeeze(sarImage)),'FaceColor','interp','LineStyle','none')
% view(2)

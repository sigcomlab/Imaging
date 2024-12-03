% Mutual Coupling Compensation for Compact MIMO Radar Article
clc
clear
%% radar properties
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 78.5e9;%77e9; % start frequency
% N = 64; % number of symbols
% N_FFT_kx = 512; % number of symbols in x-axis
% N_FFT_ky = 1024; % number of symbols in y-axis
% 
mu = 2.573e13;%70.295e12; % slope
fs = 2e6; %5e6;        % Sampling rate (sps)
Ts = 1/fs;          % Sampling period
km = mu / c;
k = 2*pi*f_0/c;

R = 1.1;
%% TX RX properties
y_T = [0 0]-102.2e-3; % 2-TX 
x_T = [0.0117 0.0195]; 
n_T = length(x_T);
x_Tc = mean (x_T);
y_Tc = mean(y_T);

y_R = [0 0 0 0]-102.2e-3;
x_R = [0  0.0019 0.0039 0.0058]; 
n_R = length(x_R);
x_Rc = mean(x_R); 
y_Rc = mean(y_R);

plot(x_T, y_T  ,'O'), hold on, grid on
plot(x_R, y_R ,'*')
plot(x_Tc, y_Tc,'+')
plot(x_Rc, y_Rc,'.')

x_c = (x_Tc + x_Rc) /2;
y_c = (y_Tc + y_Rc) /2;

C_T = eye(n_T);
C_R = eye(n_R);
theta = [-40:1:40];
K = length(theta);
khi = zeros(n_T, n_R, K);
 % eq. 4
W = 0.05*(randn(n_T,n_R)+1i*randn(n_T,n_R)); % Noise 
for q = 1: 90 
    for i = 1: K
        A = exp (-1i* ((2*pi*f_0)/c) * (x_T' + x_R) * sind(theta(i))); 
        S = exp (-1i* ((2*pi*f_0)/c) * (2*R + x_c*sind(theta(i)) + y_c*sind(theta(i)))) * A;
        X = C_T * S * C_R + W;
    
        khi(:,:,i) = X;
        phi(:,:,i) = C_T * S;
    end
    Khi = reshape(permute(khi,[2,1,3]),size(khi,2),[])';
    Phi = reshape(permute(phi,[2,1,3]),size(phi,2),[])';
    C_R = (Phi' * Phi)^-1 * Phi' * Khi;  
    
    for i = 1: K
        A = exp (-1i* ((2*pi*f_0)/c) * (x_T' + x_R) * sind(theta(i))); 
        S = exp (-1i* ((2*pi*f_0)/c) * (2*R + x_c*sind(theta(i)) + y_c*sind(theta(i)))) * A;
        X_p = C_T * S * C_R + W;
    
        khi_p(:,:,i) = X_p;
        gamma(:,:,i) = S * C_R;
    end
    Khi_p = reshape(permute(khi_p,[1,2,3]),size(khi_p,1),[]);
    Gamma = reshape(permute(gamma,[1,2,3]),size(gamma,1),[]);
    C_T = Khi_p * Gamma' * (Gamma * Gamma')^-1;

end
%% 
abs(C_T)/max(max(abs(C_T)))
angle(C_T)

abs(C_R)/max(max(abs(C_R)))
angle(C_R)



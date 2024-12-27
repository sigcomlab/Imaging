% Mutual Coupling Compensation for Compact MIMO Radar Article
clc
clear
%% Radar properties
c = 299792458; % physconst('lightspeed'); in m/s
f_0 = 78.5e9;%77e9; % start frequency
mu = 2.573e13;%70.295e12; % slope
fs = 2e6; %5e6;     % Sampling rate (sps)
Ts = 1/fs;          % Sampling period
km = mu / c;
k = 2*pi*f_0/c;
R = 1.1; % The distance between target and radar
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
theta = (-40:1:40);
% load("per_elahe_RP.mat")
% theta=squeeze(theta);
K = length(theta);
khi = zeros(n_T, n_R, K);

W = 0.01*(randn(n_T,n_R) + 1i*randn(n_T,n_R)); % Noise 

% The coupling matrixes in the article
C_Tme = [0.99*exp(1i*(0.1)) 0.0263*exp(1i*(-0.14)); 0.0199*exp(1i*(-0.01)) 1*exp(1i*(0))];  % eq. 21
C_Rme = [0.9219*exp(1i*(-0.14)) 0.0852*exp(1i*(-1.78)) 0.0781*exp(1i*(1.93)) 0.0814*exp(1i*(2.52))
         0.0712*exp(1i*(3.08)) 0.8766*exp(1i*(0.12)) 0.1211*exp(1i*(-1.74)) 0.0789*exp(1i*(-0.41))
         0.0784*exp(1i*(-0.71)) 0.1080*exp(1i*(-1.85)) 0.8107*exp(1i*(0.1)) 0.1485*exp(1i*(-3.01))
         0.0640*exp(1i*(1.73)) 0.0498*exp(1i*(1.05)) 0.1048*exp(1i*(-1.33)) 1*exp(1i*(-0))];  % eq. 22
% C_Rme = C_Rme*5 - eye(4);
% X = permute(X,[2,3,1]);
for q = 1: 90 
    for i = 1: K
        A = exp (-1i* ((2*pi*f_0)/c) * (x_T' + x_R) * sind(theta(i))); 
        S = exp (-1i* ((2*pi*f_0)/c) * (2*R + x_c*sind(theta(i)) + y_c*sind(theta(i)))) * A;
        X = C_Tme * S * C_Rme + W;
        khi(:,:,i) = X;
        phi(:,:,i) = C_T * S;
        gamma(:,:,i) = S * C_R;

        %% Plotting other plots in the article
        % SU = unwrap(angle(reshape(S.', [1,8])));
        % XU = unwrap(angle(reshape(X.', [1,8])));
        % muS(i) = SU(8)-SU(1); % Slpoe
        % muX(i) = XU(8)-XU(1);
        % if (theta(i) == 20 && q==90) % Theta = 20 was plotted in the article
        %     figure()
        %     ang_scale = asin(linspace(-1,1,64));
        %     plot(ang_scale,log10(abs(fftshift(fft(reshape(S.', [1,8]), 64)))))
        % end

        % if q == 90
        %     gamma_start = sin(theta(i)/180*3.14)
        %     gamma_first_lobe = gamma_start - 2/8;
        %     if gamma_first_lobe < -1
        %         gamma_first_lobe = gamma_first_lobe + 2;
        %     end
        %     gamma_first_lobe
        % 
        %     theta_first_lobe = asin(gamma_first_lobe)/3.14*180;
        %     [difference, theta_idx] = min(theta-theta_first_lobe);
        % 
        %     beampattern = 10*log10(abs(fftshift(fft(reshape(S.', [1,8]), 128))));
        % 
        %     sll(i) = beampattern(theta_idx)-max(beampattern)
            % figure()
            %plot(beampattern-max(beampattern))
        % end
        %%
    end
    
    Khi = reshape(permute(khi,[2,1,3]),size(khi,2),[]).';
    Phi = reshape(permute(phi,[2,1,3]),size(phi,2),[]).';
    C_R = (Phi' * Phi)^-1 * Phi' * Khi;   % eq. 17

    Khi_p = reshape(permute(khi,[1,2,3]),size(khi,1),[]);
    Gamma = reshape(permute(gamma,[1,2,3]),size(gamma,1),[]);
    C_T = Khi_p * Gamma' * (Gamma * Gamma')^-1;  % eq. 19

end
% figure(), hold on;
% plot(theta, muX - muS);
% plot(theta, sll)
% plot(muX);

%% 
abs_CT = abs(C_T)/max(max(abs(C_T)))
ang_CT = angle(C_T)

abs_CR = abs(C_R)/max(max(abs(C_R)))
ang_CR = angle(C_R)




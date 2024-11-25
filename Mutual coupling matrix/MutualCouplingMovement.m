tic
clc
clear all
clear figure
Np=6; % position number on the rail
NRx=8; % Rx number
T=56e-6; % The width of the chirp
fs=5e6; % Sampling frequency
Ns=256; % Sample number
NT=6;  % Real target number
Nx=10*NT;  % hypothetical target number
I = eye(Np*Ns); % [Np*Ns , Np*Ns]
mu=87e12; % Slope of chirp
k=1; % The variable in constraint equation (34b)
C0=eye(NRx); % Primary matrix for coupling matrix in (33b)
f0=77e9; % center frequency
c = physconst('Lightspeed'); % speed velocity
lambda=c/f0; % Wavelength
% rad = 0.040; AT = mean(mean(rcssphere(rad,c,f0)));
radar_step =- 0.04; % radar step length
AT=0.04; % Rdar RCS 
%% TX RX target location
% TX location defenition
xt=0;
yt=0;
zt=-28e-3;
% First RX1 location x_target_position
xr1= -2e-3;
yr=0; 
zr=28e-3;
%% Hypothetical_target_position between 1.5 to 2.5 m
x_target_position = linspace(1*lambda,-1*NRx*Np*lambda,Nx); % x-target location 
y_target_position = 1.5+rand(1, Nx);
z_target_position = zeros(1, Nx);
scatter3(x_target_position,y_target_position,z_target_position,"filled"),  hold on, grid on, 
%% Real target location
real_point=randperm(length(x_target_position), NT);
x_target_position_real = x_target_position(real_point)%+1e-2*x_target_position(real_point);
y_target_position_real = y_target_position(real_point)%+1e-2* y_target_position(real_point);
z_target_position_real = z_target_position(real_point)%+1e-2*z_target_position(real_point);
scatter3(x_target_position_real,y_target_position_real,z_target_position_real,"filled")
hold on, 
%% H matrix calculation in equation 17
[tau,aTX,xt_plot,xr_plot] = delay_attenuation_calculation(xt,yt,zt,xr1,yr,zr,x_target_position,y_target_position,z_target_position,NRx,lambda,Np,radar_step);
H = H_calculation (Nx ,Np, Ns, NRx,f0,tau,aTX,mu,T);
%% plot TX and RX on the rail
scatter3(xt_plot,yt*ones(1,length(xt_plot)),zt*ones(1,length(xt_plot)));
scatter3(xr_plot,yr*ones(length(xr_plot),1),zr*ones(length(xr_plot),1)); 
xlabel('x'),ylabel('y'),zlabel('z')
legend('Hypothetical-target-position','Real-target-position','TX-position-on-the-rail','RX-position-on-the-rail')
%%  Simulation
load 'C_article' 
%% cration of measured data y in equation (23)
Cy=C_article;  % coupling matrix reported by article
kron_I_Cy=kron(I, Cy);
[tau_m,aTX_m] = delay_attenuation_calculation(xt,yt,zt,xr1,yr,zr,x_target_position_real,y_target_position_real,z_target_position_real,NRx,lambda,Np,radar_step);
H_m = H_calculation (NT ,Np, Ns, NRx,f0,tau_m,aTX_m,mu,T);
y = kron_I_Cy*H_m*AT * ones(NT,1); 
figure % plotting one sample of created signal
plot([1:Ns]*(c/(2*mu))*fs/Ns,abs(y(1:8:Ns*8)))
grid on
Ey=norm(y,2)^2;
%%  Optimization problem in Algorithm 1 (34a)
for iteration=1:10
    cvx_begin
    variable x(Nx,1) 
    kron_I_C=kron(I, C0);
    minimize(norm(kron_I_C*H*x - y ,2))
       subject to
            real(x) >= 0;
           (norm(x,1)) <= k*NT*AT; % constraint in (34b)
    cvx_end
%% estimate_position (xi)
[sorted_vec, sorted_idx] = sort(x,'descend');
p_estimate_position_x = x_target_position(sorted_idx(1:Nx));
p_estimate_position_y = y_target_position(sorted_idx(1:Nx));
p_estimate_position_z = z_target_position(sorted_idx(1:Nx));
[tau_m,aTX_m] = delay_attenuation_calculation(xt,yt,zt,xr1,yr,zr,p_estimate_position_x,p_estimate_position_y,p_estimate_position_z,NRx,lambda,Np,radar_step);
H_m1 = H_calculation (NT ,Np, Ns, NRx,f0,tau_m,aTX_m,mu,T);
x = sorted_vec (1:NT);
%% Problem IN Algorithm 1 (34c)
    cvx_begin
    variable C_prime(NRx,NRx)
      subject to
          for i=1: length(C_prime)
            for j=1: length(C_prime)
                     C_prime (i,j) == C_prime (j,i);
                     C_prime (i,i)==1; % equality of diagonal element of mutual coupling matrix
                     C_prime (i,j)>0;
                     C_prime (i,j)<=1;
             end
          end
     kron_I_C_prime=kron(I, C_prime);
     minimize(norm(kron_I_C_prime*H_m1*x - y ,2)) 
    cvx_end
kron_I_C_prime = kron(I, C_prime);
Eh = norm(kron_I_C_prime*H_m1*x ,2)^2;
C = sqrt (Ey/Eh) * C_prime,
C0=C;
CC(:,:,iteration)=C;
XX(:,iteration)=x;
end
toc
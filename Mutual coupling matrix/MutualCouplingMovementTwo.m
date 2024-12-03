% Antenna Array Calibration Using a Sparse Scene Article
% 30/05/2023 General NP NT NRx NS
% article primary data for creating y
% it works only by deviation 1e-5*()
tic
clc
clear all
clear figure
Np=6; % position number on the rail
NRx=8; % Rx number
% B=4e9;
T=56e-6;
fs=5e6;
Ns=256; % Sample number
NT=6; 
Nx=10*NT; 
I = eye(Np*Ns); % [Np*NRx*Ns , Np*NRx*Ns]
mu=87e12; %slope
k=1;
C0=eye(NRx);
f0=77e9; c = physconst('Lightspeed'); lambda=c/f0;
% rad = 0.040; AT = mean(mean(rcssphere(rad,c,f0)));
radar_step =- 0.04;
AT=4;
%% TX RX target location
% TX location
xt=0;
yt=0;
zt=-28e-3;
% First RX1 location x_target_position
xr1= -2e-3;
yr=0; 
zr=28e-3;
%% Hypothetical_target_position
x_target_position = linspace(1*lambda,-1*NRx*Np*lambda,Nx); % x-target location 
y_target_position = 1.5+rand(1, Nx);
z_target_position = zeros(1, Nx);
e = Nx/2;
x_target_position (e)= -0.12; % x-target location 
y_target_position (e)= 0.908584;
z_target_position (e)= 0;
scatter3(x_target_position,y_target_position,z_target_position),  hold on, grid on
scatter3(x_target_position(Nx/2),y_target_position(Nx/2),z_target_position(Nx/2),"filled"),  hold on, grid on
%% One Real target location
real_point=randperm(length(x_target_position), NT);
% real_point=(e);
x_target_position_real = x_target_position(real_point)%+1e-5*x_target_position(real_point);
y_target_position_real = y_target_position(real_point)%+1e-5* y_target_position(real_point);
z_target_position_real = z_target_position(real_point)%+1e-5*z_target_position(real_point);
scatter3(x_target_position_real,y_target_position_real,z_target_position_real,"filled")
hold on
%% radar measurment
[tau,aTX,xt_plot,xr_plot] = delay_attenuation_calculation(xt,yt,zt,xr1,yr,zr,x_target_position,y_target_position,z_target_position,NRx,lambda,Np,radar_step);
H = H_calculation (Nx ,Np, Ns, NRx,f0,tau,aTX,mu,T);
%% plot TX and RX on the rail
scatter3(xt_plot,yt*ones(1,length(xt_plot)),zt*ones(1,length(xt_plot)));
scatter3(xr_plot,yr*ones(length(xr_plot),1),zr*ones(length(xr_plot),1));
xlabel('x'),ylabel('y'),zlabel('z')
%%  Simulation
load 'C_article' 
% C_article=(eye(NRx)); 
%% applying exp (j*pi)
% for i = 1: NRx
%     for j = 1: NRx
%        C_article (i,j) = C_article (i,j) .* exp (1i*(i+j)*pi);
%     end
% end
%%
Cy=C_article;  
kron_I_Cy=kron(I, Cy);
[tau_m,aTX_m] = delay_attenuation_calculation(xt,yt,zt,xr1,yr,zr,x_target_position_real,y_target_position_real,z_target_position_real,NRx,lambda,Np,radar_step);
H_m = H_calculation (NT ,Np, Ns, NRx,f0,tau_m,aTX_m,mu,T);
y_article = kron_I_Cy*H_m*AT * ones(NT,1); 
y = y_article;
yme=y;
figure
plot([1:Ns]*(c/(2*mu))*fs/Ns,abs(yme(1:8:Ns*8)))
grid on
%%
%%  data=data(1:40:318,:,8,9:16); 03-04-2023_17-02-01__one_sphere_88_cm.npz
% load data.mat  % size[8,512,1,8]
% A=squeeze(data);
% figure
% % plot(abs(A(1,:,1)))
% hold on
% for ks=1:8
%     for kk=1:8
%         A_clean= A(ks,:,kk);
%         for i = 1: length(A_clean)
%             if i<22 || i>25
% %             if abs(A_clean(i))<3e4 || i>300
% 
%                   A_clean(i)=0;
%             end
%         end
%         A(ks,:,kk)=A_clean;
% %         plot(abs(A(ks,:,kk)))
%     end
% end
% % plot(abs(A(6,:,4)))
% B = permute(A,[3 2 1]);
% c1 = squeeze(num2cell(B,[1]));
% c2 = c1(:);
% % y =18*cat(1, c2{:})/6;
% plot([1:Ns]*(c/(2*mu)*fs)/Ns,abs(y(1:8:Ns*8)))
% U=[1:Ns]*(c/(2*mu)*fs)/Ns;
 Ey=norm(y,2)^2;
%% Second Problem
for iteration=1:6
    cvx_begin
    variable x(Nx,1) 
    kron_I_C=kron(I, C0);
    minimize(norm(kron_I_C*H*x - y ,2))
       subject to
            real(x) >= 0;
           (norm(x,1)) <= k*NT*AT;
    cvx_end
%% estimate_position (xi)
[sorted_vec, sorted_idx] = sort(x,'descend')
sum(x)
p_estimate_position_x = x_target_position(sorted_idx(1:Nx));
p_estimate_position_y = y_target_position(sorted_idx(1:Nx));
p_estimate_position_z = z_target_position(sorted_idx(1:Nx));
[tau_m,aTX_m] = delay_attenuation_calculation(xt,yt,zt,xr1,yr,zr,p_estimate_position_x,p_estimate_position_y,p_estimate_position_z,NRx,lambda,Np,radar_step);
H_m1 = H_calculation (NT ,Np, Ns, NRx,f0,tau_m,aTX_m,mu,T);
x = sorted_vec (1:NT);
%% First Problem
    cvx_begin
    variable C_prime(NRx,NRx)
      subject to
          for j = NRx:-1:1
                for i = 1:j-2
                      C_prime(i,i+(NRx-(j-1))) == C_prime(i+1,i+(NRx-(j-2)));
                end
          end
          for i=1: length(C_prime)
                for j=1: length(C_prime)
%                     if mod(i+j, 2) == 1
%                          C_prime((i-1)*NRx + j) <= -1;  % Set negative sign constraint
%                     else
%                           C_prime((i-1)*NRx + j) >= 1;   % Set positive sign constraint
%                     end
%                    C_prime (i,j) == C_prime (i,j)* exp (1i*(i+j)*pi);
                     C_prime (i,j) == C_prime (j,i);
                     C_prime (i,i)==1;
                     C_prime (i,j)>0;
                     C_prime (i,j)<=1;
                end
          end
     kron_I_C_prime=kron(I, C_prime);
     minimize(norm(kron_I_C_prime*H_m1*x - y ,2))
    cvx_end
kron_I_C_prime = kron(I, C_prime);
C_prime
Eh = norm(kron_I_C_prime*H_m1*x ,2)^2;
C = sqrt (Ey/Eh) * C_prime,
C0=C;
CC(:,:,iteration)=C;
XX(:,iteration)=x
end
toc

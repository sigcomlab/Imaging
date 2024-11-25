function H = H_calculation(Nx ,Np,Ns,NRx,f0,tau,aTX,mu,T)
% H1=zeros(NRx,Ns);  % Np.NRx.Ns
% H2=zeros(NRx*Ns,Np);  % Np.NRx.Ns
% H=zeros(NRx*Ns*Np,Nx);  % Np.NRx.Ns
% for nx=1:Nx
%     for np=1:Np
%         for ns=1:Ns
%             for nrx=1:NRx
%                  H1(nrx,ns)=exp(1i*2*pi*tau(np,nrx)*(mu*ns*T/Ns + f0));
%             end
%         end
%         H2(:,np)=reshape(H1,[Ns*NRx,1]);
%     end
% H(:,nx)=reshape(H2,[Np*Ns*NRx,1])/sqrt(Ns);
% end
%% delay-tau calculation
% % TX location
% xt=0;
% yt=0;
% 
% % First RX1 location
% xr1=lambda;
% yr=0; 
% 
% x_target_position = linspace(0,5*lambda/2,Nx); % x-target location 
% y_target_position = 1 * ones(1,Nx); % y-target location 
% %% RT , RR Calculation
% for i =1:length(x_target_position)
%         RT(i)=sqrt((x_target_position(i)-xt)^2+(y_target_position(i)-yt)^2); 
%     for j=0:NRx-1
%             RR(j+1,i)=sqrt((x_target_position(i)-(xr1+j*lambda/2))^2+(y_target_position(i)-yr)^2); 
%             tau(j+1,i)=(RT(i) + RR(j+1))/c;
%             aTX(j+1,i) = 1/(RT(i)*RR(j+1));
%     end    
% end
% tau = tau';
% aTX = aTX';
%% W-Calculation
for i=1:Ns
    for j=1:Ns
        W(i,j) = exp(-1i*2*pi/Ns)^((i-1)*(j-1))/sqrt(Ns);
    end
end
%% 
S_id_np_nRX = zeros(Ns,Nx);
for np=1:Np
for nrx=1:NRx
    for nx=1:Nx
        for ns=1:Ns
                
                S_id_np_nRX(ns,nx)=aTX(nx,nrx,np)*exp(1i*2*pi*tau(nx,nrx,np)*(mu*(ns-1)*T/Ns + f0));
        end
    end
        S_DFT_id_np_nRX = W * S_id_np_nRX; 
        A (:,:,nrx) = S_DFT_id_np_nRX ;
end


for nx =1:Nx
    for ns=1:Ns
        for nrx =1 : NRx
            B(nrx,ns) = A(ns,nx,nrx); 
        end
    end
     H_prime_np_ns (:,nx)  = reshape(B,[NRx*Ns,1]);
end
H_np (:,:,np)  = H_prime_np_ns;
end
H = squeeze(num2cell(H_np,[1,2])); 
H=vertcat(H{:});
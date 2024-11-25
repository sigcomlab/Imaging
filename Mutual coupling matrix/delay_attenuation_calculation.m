function [tau,aTX,xt_plot,xr_plot] = delay_attenuation_calculation(xt,yt,zt,xr1,yr,zr,x_target_position,y_target_position,z_target_position,NRx,lambda,Np,radar_step)
c = physconst('Lightspeed');
for k = 0:Np-1
    for i =1:length(x_target_position)
               RT(k+1,i)=sqrt((x_target_position(i)-(xt+k*radar_step))^2+(y_target_position(i)-yt)^2)+(z_target_position(i)-zt)^2; 
               xt_plot(k+1) = xt+k*radar_step;
        for j=0:NRx-1
               RR(j+1,i,k+1)=sqrt((x_target_position(i)-(xr1+j*lambda/2+k*radar_step))^2+(y_target_position(i)-yr)^2)+(z_target_position(i)-zr)^2; 
               tau(j+1,i,k+1)=(RT(k+1,i) + RR(j+1,i,k+1))/c;
               aTX(j+1,i,k+1) = 1/(RT(k+1,i)*RR(j+1,i,k+1));
               xr_plot(j+1,k+1) = xr1+j*lambda/2+k*radar_step;
        end 
    end
end

xr_plot = reshape(xr_plot,[],1);
% for i =1:length(x_target_position)
%         RT(i)=sqrt((x_target_position(i)-xt)^2+(y_target_position(i)-yt)^2); 
%     for j=0:NRx-1
%             RR(j+1,i)=sqrt((x_target_position(i)-(xr1+j*lambda/2))^2+(y_target_position(i)-yr)^2); 
%             tua(j+1,i)=(RT(i) + RR(j+1))/c;
%             aTX(j+1,i) = 1/(RT(i)*RR(j+1));
%     end    
% end
tau = permute (tau,[2,1,3]);
aTX = permute (aTX,[2,1,3]);
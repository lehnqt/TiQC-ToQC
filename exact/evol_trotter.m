function [I,ovl]=evol_trotter(H0,Hc,I0,time_grid,c0,c)
cert_num=length(H0);
ctrl_num=length(Hc);
bin_num=length(time_grid)-1;
c0=reshape(c0,[bin_num,cert_num]);
c=reshape(c,[bin_num,ctrl_num]);
I=I0;
dt=time_grid(2)-time_grid(1);
if nargout<2
for j=1:bin_num
    %if mod(j-1,100)==0
        %fprintf('step %d\n',j);
    %end
    H_int=sparse(0);
    for k=1:cert_num
     H_int=H_int+H0(k).op*H0(k).ft(time_grid(j)+dt/2)*c0(j,k);
    end
    H_tot=sparse(0);
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*c(j,k);
    end
    I=expm(-1i*dt*H_tot)*expm(-1i*dt*H_int)*I*expm(1i*dt*H_int)*expm(1i*dt*H_tot);
end
else
    n=log2(length(I0));
    ovl=zeros(bin_num,1);
    normI=n*2^n;
for j=1:bin_num
    % if mod(j-1,100)==0
    %     fprintf('step %d\n',j);
    % end
    H_int=sparse(0);
    for k=1:cert_num
     H_int=H_int+H0(k).op*H0(k).ft(time_grid(j)+dt/2)*c0(j,k);
    end
    H_tot=sparse(0);
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*c(j,k);
    end
    I=expm(-1i*dt*H_tot)*expm(-1i*dt*H_int)*I*expm(1i*dt*H_int)*expm(1i*dt*H_tot);
    ovl(j)=trace(I0'*I)/normI;
end
end
end
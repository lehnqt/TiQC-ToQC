function I=evol(H0,Hc,I0,time_grid,c0,c)
cert_num=length(H0);
ctrl_num=length(Hc);
bin_num=length(time_grid)-1;
c0=reshape(c0,[bin_num,cert_num]);
c=reshape(c,[bin_num,ctrl_num]);
I=I0;
for j=1:bin_num
    dt=time_grid(j+1)-time_grid(j);
    H_tot=sparse(0);
    for k=1:cert_num
        H_tot=H_tot+H0(k).op*H0(k).ft(time_grid(j)+dt/2)*c0(j,k);
    end
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*c(j,k);
    end
    I=expm(-1i*dt*H_tot)*I;
    I=expm(-1i*dt*H_tot)*I';
    I=I';
end
end
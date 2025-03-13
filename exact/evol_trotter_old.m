function I=evol_trotter(H0,Hc,I0,time_grid,c,iso)
cert_num=length(H0);
ctrl_num=length(Hc);
bin_num=length(time_grid)-1;
c=reshape(c,[bin_num,ctrl_num]);
I=I0;
H_int=sparse(0);
for k=1:cert_num
     H_int=H_int+H0(k).op;
end
dt=time_grid(2)-time_grid(1);
if iso==1
 I=expm(1i*(dt/2)*H_int)*I*expm(-1i*(dt/2)*H_int);   
end
for j=1:bin_num
    H_tot=sparse(0);
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*c(j,k);
    end
    I=expm(-1i*dt*H_tot)*expm(-1i*dt*H_int)*I*expm(1i*dt*H_int)*expm(1i*dt*H_tot);
end
if iso==1
 I=expm(-1i*(dt/2)*H_int)*I*expm(1i*(dt/2)*H_int);   
end
end
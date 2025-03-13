function [iF,iG,overlap]=infid_krylov(H0,Hc,I0,Itg,time_grid,c,numK)
cert_num=length(H0);
ctrl_num=length(Hc);
bin_num=length(time_grid)-1;
c=reshape(c,[bin_num,ctrl_num]);
I=I0;
if nargout<2
for j=1:bin_num
    dt=time_grid(j+1)-time_grid(j);
    H_tot=sparse(0);
    for k=1:cert_num
        H_tot=H_tot+H0(k).op*H0(k).ft(time_grid(j)+dt/2);
    end
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*c(j,k);
    end
    I=expv_krylov(dt,H_tot,I,numK);
    I=expv_krylov(dt,H_tot,I',numK);
    I=I';
end
    Ovl=trace(Itg'*I);
    iF=1-real(Ovl);
else
iG=zeros(bin_num,ctrl_num);
overlap=zeros(bin_num,1);
D=cell(bin_num,ctrl_num);
Ifw=cell(bin_num+1,1);
Ibw=cell(bin_num+1,1);
Ifw{1}=I0;
Ibw{1}=Itg;
for j=1:bin_num
    %forward propagation
    dt=time_grid(j+1)-time_grid(j);
    H_tot=sparse(0);
    for k=1:cert_num
        H_tot=H_tot+H0(k).op*H0(k).ft(time_grid(j)+dt/2);
    end
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*c(j,k);
    end
    Itemp=expv_krylov(dt,H_tot,Ifw{j},numK);
    Itemp=expv_krylov(dt,H_tot,Itemp',numK);
    Ifw{j+1}=Itemp';
    for k=1:ctrl_num
        C1=Hc(k).op*Hc(k).ft(time_grid(j)+dt/2);
        C2=H_tot*C1-C1*H_tot;
        D{j,k}=-1i*dt*C1+(dt^2/2)*C2;
%         C3=(H0{i}+H_int)*C2-C2*(H0{i}+H_int);
%         D{k,j}=-1i*Delta_t*C1+(Delta_t^2/2)*C2+1i*(Delta_t^3/6)*C3;
    end
    %backward propagation
    dt=time_grid(bin_num+2-j)-time_grid(bin_num+1-j);
    H_tot=sparse(0);
    for k=1:cert_num
        H_tot=H_tot+H0(k).op*H0(k).ft(time_grid(bin_num+1-j)+dt/2);
    end
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(bin_num+1-j)+dt/2)*c(bin_num+1-j,k);
    end
   Itemp=expv_krylov(-dt,H_tot,Ibw{j},numK);
   Itemp=expv_krylov(-dt,H_tot,Itemp',numK);
   Ibw{j+1}=Itemp';
end
    Ovl=trace(Ibw{1}'*Ifw{bin_num+1});
    iF=1-real(Ovl);
for j=1:bin_num
    overlap(j)=trace((Ibw{bin_num-j+2})'*Ifw{j});
    for k=1:ctrl_num
        iG(j,k)=-real(trace(Ibw{bin_num-j+2}'*(D{j,k}*Ifw{j}+Ifw{j}*(D{j,k})')));
    end
end
iG=iG(:);
end
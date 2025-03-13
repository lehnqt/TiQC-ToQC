function [iFlist,iGlist,overlap]=robust_infid_exact_sgd(H0,Hc,Hu,err,sampl_num,I0,Itg,time_grid,c,ismean)
cert_num=length(H0);
ctrl_num=length(Hc);
unc_num=length(Hu);
u=err*(rand(sampl_num,unc_num)-1/2);
bin_num=length(time_grid)-1;
c=reshape(c,[bin_num,ctrl_num]);
dt=time_grid(2)-time_grid(1);
if nargout<2
    iFlist=zeros(sampl_num,1);
    for js=1:sampl_num
    ur=u(js,:);
    H_const=sparse(0);
    for k=1:cert_num
        H_const=H_const+H0(k).op;
    end
    for k=1:unc_num
        H_const=H_const+ur(k)*Hu(k).op;
    end
    I=I0;
for j=1:bin_num
    H_tot=H_const;
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*c(j,k);
    end
    I=expm(-1i*dt*H_tot)*I;
    I=expm(-1i*dt*H_tot)*I';
    I=I';
end
    Ovl=trace(Itg'*I);
    iFlist(js)=1-real(Ovl);
    end
    if ismean==1
        iFlist=mean(iFlist);
    end
else
    iGlist=zeros(bin_num*ctrl_num,sampl_num);
    iFlist=zeros(sampl_num,1);
    for js=1:sampl_num
iG=zeros(bin_num,ctrl_num);
overlap=zeros(bin_num,1);
D=cell(bin_num,ctrl_num);
Ifw=cell(bin_num+1,1);
Ibw=cell(bin_num+1,1);
Ifw{1}=I0;
Ibw{1}=Itg;
  ur=u(js,:);
H_const=sparse(0);
    for k=1:cert_num
        H_const=H_const+H0(k).op;
    end
    for k=1:unc_num
        H_const=H_const+ur(k)*Hu(k).op;
    end
for j=1:bin_num
    %forward propagation
    H_tot=H_const;
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*c(j,k);
    end
    Itemp=expm(-1i*dt*H_tot)*Ifw{j};
    Itemp=expm(-1i*dt*H_tot)*Itemp';
    Ifw{j+1}=Itemp';
    for k=1:ctrl_num
        C1=Hc(k).op;
        C2=H_tot*C1-C1*H_tot;
        D{j,k}=-1i*dt*C1+(dt^2/2)*C2;
%         C3=(H0{i}+H_int)*C2-C2*(H0{i}+H_int);
%         D{k,j}=-1i*Delta_t*C1+(Delta_t^2/2)*C2+1i*(Delta_t^3/6)*C3;
    end
    %backward propagation
    dt=time_grid(bin_num+2-j)-time_grid(bin_num+1-j);
    H_tot=H_const;
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*c(bin_num+1-j,k);
    end
   Itemp=expm(1i*dt*H_tot)*Ibw{j};
   Itemp=expm(1i*dt*H_tot)*Itemp';
   Ibw{j+1}=Itemp';
end
    Ovl=trace(Ibw{1}'*Ifw{bin_num+1});
    iFlist(js)=1-real(Ovl);
for j=1:bin_num
    overlap(j)=trace((Ibw{bin_num-j+2})'*Ifw{j});
    for k=1:ctrl_num
        iG(j,k)=-real(trace(Ibw{bin_num-j+2}'*(D{j,k}*Ifw{j}+Ifw{j}*(D{j,k})')));
    end
end
iGlist(:,js)=iG(:);
    end
    iFlist=mean(iFlist);
    iGlist=mean(iGlist,2);
end
end
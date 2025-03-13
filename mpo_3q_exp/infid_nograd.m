function [iF,maxD]=infid_nograd(H0,Hc,c,T,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso)
n=length(mpo0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
d=size(mpo0{1},2);
dt=T/(nbin*nt);
maxD=1;
if iso==1
%initial half time evolution
g20=cell(1,n-1);
for j=1:n-1
        h=H0{j};
        gate=expm(1i*dt*h/2);%half time, backward
        gate=reshape(gate,[d,d,d,d]);
        g20{j}=gate;
end
%apply odd terms
for j=1:2:n-1
    [mpo0{j},mpo0{j+1}]=gate_2q_LR(mpo0{j},mpo0{j+1},g20{j},sv_min,D);
end
 %apply even terms
for j=2:2:n-1
    [mpo0{j},mpo0{j+1}]=gate_2q_LR(mpo0{j},mpo0{j+1},g20{j},sv_min,D);
end
%backward
%apply even terms
for j=2:2:n-1
    [mpotg{j},mpotg{j+1}]=gate_2q_LR(mpotg{j},mpotg{j+1},g20{j},sv_min,D);
end
%apply odd terms
for j=1:2:n-1
    [mpotg{j},mpotg{j+1}]=gate_2q_LR(mpotg{j},mpotg{j+1},g20{j},sv_min,D);
end
if iscpr==1
      mpo0=mpo_recompress(mpo0,sv_min,Dc,nsweep);
      mpotg=mpo_recompress(mpotg,sv_min,Dc,nsweep);
end
mpo0=mpo_normalize(mpo0);
mpotg=mpo_normalize(mpotg);
end
%normal evolution
mpo=mpo0;
g2=cell(1,n-1);
for j=1:n-1
        h=H0{j};
        gate=expm(-1i*dt*h);
        gate=reshape(gate,[d,d,d,d]);
        g2{j}=gate;
end
for k=1:nbin
     %gate construction
   g1=cell(1,n);
    for j=1:n
       h=zeros(d);
       for jc=1:nc
            for js=1:length(Hc(jc).sys)
                if Hc(jc).sys(js)==j
        h=h+c(k,jc)*Hc(jc).op{js};
                end
            end
       end
       gate=expm(-1i*dt*h);
       g1{j}=gate;
    end 
%forward propagation
for jt=1:nt
 %apply odd terms
for j=1:2:n-1
    [mpo{j},mpo{j+1}]=gate_2q_LR(mpo{j},mpo{j+1},g2{j},sv_min,D);
end
 %apply even terms
for j=2:2:n-1
    [mpo{j},mpo{j+1}]=gate_2q_LR(mpo{j},mpo{j+1},g2{j},sv_min,D);
end
%apply 1q terms
for j=1:n
    [mpo{j}]=gate_1q_LR(mpo{j},g1{j});
end
end
if mod(k-1,midstep)==0||k==nbin
    if iscpr==1
        mpo=mpo_recompress(mpo,sv_min,Dc,nsweep);
    end
mpo=mpo_normalize(mpo);
end
for j=1:n
maxD=max(maxD,max(size(mpo{j})));
end
end
ovl=overlap(mpotg,mpo);
iF=1-real(ovl);   
end

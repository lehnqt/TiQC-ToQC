function [iFlist,Dlist]=mps_infid_nograd_robust(J,H0,Hc,c,T,mps0,mpstg,sv_min,D,midstep,nt,iscpr,iso,ismean)
n=length(mps0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
d=size(mps0{1},2);
dt=T/(nbin*nt);
nsampl=size(J,1);
iFlist=zeros(nsampl,1);
Dlist=zeros(nsampl,1);
for jsampl=1:nsampl
maxD=1;
if iso==1
    %initial half time evolution
g20=cell(1,n-1);
for j=1:n-1
        h=J(jsampl,j)*H0{j};
        gate=expm(1i*dt*h/2);%half time, backward
        gate=reshape(gate,[d,d,d,d]);
        g20{j}=gate;
end
%apply odd terms
for j=1:2:n-1
    [mps0{j},mps0{j+1}]=mps_gate_2q(mps0{j},mps0{j+1},g20{j},sv_min,D);
end
 %apply even terms
for j=2:2:n-1
    [mps0{j},mps0{j+1}]=mps_gate_2q(mps0{j},mps0{j+1},g20{j},sv_min,D);
end
%backward
%apply even terms
for j=2:2:n-1
    [mpstg{j},mpstg{j+1}]=mps_gate_2q(mpstg{j},mpstg{j+1},g20{j},sv_min,D);
end
%apply odd terms
for j=1:2:n-1
    [mpstg{j},mpstg{j+1}]=mps_gate_2q(mpstg{j},mpstg{j+1},g20{j},sv_min,D);
end
if iscpr==1
      mps0=mps_recompress(mps0,sv_min,D,2);
      mpstg=mps_recompress(mpstg,sv_min,D,2);
end
mps0=mps_normalize(mps0);
mpstg=mps_normalize(mpstg);
end
mps=mps0;
g2=cell(1,n-1);
for j=1:n-1
        h=J(jsampl,j)*H0{j};
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
    [mps{j},mps{j+1}]=mps_gate_2q(mps{j},mps{j+1},g2{j},sv_min,D);
end
 %apply even terms
for j=2:2:n-1
    [mps{j},mps{j+1}]=mps_gate_2q(mps{j},mps{j+1},g2{j},sv_min,D);
end
%apply 1q terms
for j=1:n
    [mps{j}]=mps_gate_1q(mps{j},g1{j});
end
end
if mod(k-1,midstep)==0||k==nbin
    if iscpr==1
        mps=mps_recompress(mps,sv_min,D,2);
    end
mps=mps_normalize(mps);
end
for j=1:n
maxD=max(maxD,max(size(mps{j})));
end
end
ovl=mps_overlap(mpstg,mps);
iF=1-abs(ovl)^2;   
iFlist(jsampl)=iF;
Dlist(jsampl)=maxD;
end
if ismean==1
    iFlist=mean(iFlist);
    Dlist=mean(Dlist);
end
end

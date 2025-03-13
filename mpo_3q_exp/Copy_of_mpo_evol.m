function [mpo,Dlist]=mpo_evol(H0,Hc,c0,c,T,mpo0,sv_min,D,midstep,nt,iscpr)
n=length(mpo0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
m0=size(H0,2);
c0=reshape(c0,[nbin,m0]);
d=size(mpo0{1},2);
dt=T/(nbin*nt);
mpo=mpo0;
Dlist=zeros(ceil(nbin/midstep)+2,1);
jD=0;
g2=cell(1,n-1);
for k=1:nbin
     %gate construction
   for j=1:n-1
       h=0;
       for l=1:m0
        h=h+H0{j,l}*c0(k,l);
       end
        gate=expm(-1i*dt*h);
        gate=reshape(gate,[d,d,d,d]);
        g2{j}=gate;
   end
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
mpo=mpo_recompress(mpo,sv_min,D,2);
    end
mpo=mpo_normalize(mpo);
jD=jD+1;
maxD=1;
for j=1:n
maxD=max(maxD,max(size(mpo{j})));
end
Dlist(jD)=maxD;
end
end   
end

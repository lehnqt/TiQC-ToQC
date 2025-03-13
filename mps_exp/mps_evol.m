function [mps,Dlist]=mps_evol(H0,Hc,c0,c,T,mps0,sv_min,D,midstep,nt)
n=length(mps0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
d=size(mps0{1},2);
dt=T/(nbin*nt);
mps=mps0;
Dlist=zeros(ceil(nbin/midstep)+2,1);
jD=0;
g2=cell(1,n-1);
for k=1:nbin
     %gate construction
   for j=1:n-1
        h=H0{j}*c0(k);
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
mps=mps_recompress(mps,sv_min,D,2);
mps=mps_normalize(mps);
jD=jD+1;
maxD=1;
for j=1:n
maxD=max(maxD,max(size(mps{j})));
end
Dlist(jD)=maxD;
end
end   
end

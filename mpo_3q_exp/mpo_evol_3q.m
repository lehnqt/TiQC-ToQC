function [mpo,Dlist]=mpo_evol_3q(H0,H1,Hc,c0,c1,c,T,mpo0,sv_min,D,midstep,nt,iscpr)
n=length(mpo0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
m0=size(H0,2);
c0=reshape(c0,[nbin,m0]);
m1=size(H1,2);
d=size(mpo0{1},2);
id=speye(d);
dt=T/(nbin*nt);
mpo=mpo0;
Dlist=zeros(ceil(nbin/midstep)+2,1);
jD=0;
g3=cell(1,n-1);
for k=1:nbin
     %gate construction
     %2q and 3q terms
  for j=1:n-2
       h=0;
       for l=1:m1
        h=h+H1{j,l}*c1(k,l);
       end
       for l=1:m0
        h=h+kron(H0{j,l}*c0(k,l),id);
       end
       if j==n-2
       for l=1:m0
          h=h+kron(H0{j+1,l}*c0(k,l),id);
       end
       end
        gate=expm(-1i*dt*h);
        gate=reshape(gate,[d,d,d,d,d,d]);
        g3{j}=gate;
  end
  %1q terms
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
%3q
%apply 1st terms
for j=1:3:n-2
    [mpo{j},mpo{j+1},mpo{j+2}]=gate_3q_LR(mpo{j},mpo{j+1},mpo{j+2},g3{j},sv_min,D);
end
 %apply 2nd terms
for j=2:3:n-2
     [mpo{j},mpo{j+1},mpo{j+2}]=gate_3q_LR(mpo{j},mpo{j+1},mpo{j+2},g3{j},sv_min,D);
end
 %apply 3rd terms
for j=3:3:n-2
    [mpo{j},mpo{j+1},mpo{j+2}]=gate_3q_LR(mpo{j},mpo{j+1},mpo{j+2},g3{j},sv_min,D);
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

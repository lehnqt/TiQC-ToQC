function [normfw,maxDfw,normbw,maxDbw,ovl_check,mpofw,mpobw]=norm_evol(H0,Hc,c,T,mpo0,mpotg,sv_min,D,midstep,nt)
n=length(mpo0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
d=size(mpo0{1},2);
dt=T/(nbin*nt);
mpofw=cell(nbin+1,n);  
mpobw=cell(nbin+1,n);
normfw=zeros(nbin,1);
normbw=zeros(nbin,1);
maxDfw=zeros(nbin,1);
maxDbw=zeros(nbin,1);
ovl_check=zeros(nbin,1);
%%%%%%%%%%%%%%%%%%%%%%%%
mpofw(1,:)=mpo0;
mpobw(1,:)=mpotg;
g2=cell(1,n-1);
for j=1:n-1
        h=H0{j};
        gate=expm(-1i*dt*h);
        gate=reshape(gate,[d,d,d,d]);
        g2{j}=gate;
end
g2bw=cell(1,n-1);
for j=1:n-1
        h=H0{j};
        gate=expm(1i*dt*h);
        gate=reshape(gate,[d,d,d,d]);
        g2bw{j}=gate;
end
for k=1:nbin
    k
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
    mpofw(k+1,:)=mpofw(k,:);
    for jt=1:nt 
     %apply odd 2q terms
    for j=1:2:n-1
        [mpofw{k+1,j},mpofw{k+1,j+1}]=gate_2q_LR(mpofw{k+1,j},mpofw{k+1,j+1},g2{j},sv_min,D);
    end
     %apply even 2q terms
    for j=2:2:n-1
        [mpofw{k+1,j},mpofw{k+1,j+1}]=gate_2q_LR(mpofw{k+1,j},mpofw{k+1,j+1},g2{j},sv_min,D);
    end
    %apply 1q terms
    for j=1:n
        [mpofw{k+1,j}]=gate_1q_LR(mpofw{k+1,j},g1{j});
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%backward propagation
%construct gates
g1=cell(1,n);
    for j=1:n
        h=zeros(d);
        for jc=1:nc
            for js=1:length(Hc(jc).sys)
                if Hc(jc).sys(js)==j
        h=h+c(nbin-k+1,jc)*Hc(jc).op{js};
                end
            end
        end
        gate=expm(1i*dt*h);
        g1{j}=gate;
    end
    mpobw(k+1,:)=mpobw(k,:);
    for jt=1:nt
    for j=1:n
        [mpobw{k+1,j}]=gate_1q_LR(mpobw{k+1,j},g1{j});
    end
    %apply even 2q terms
    for j=2:2:n-1
        [mpobw{k+1,j},mpobw{k+1,j+1}]=gate_2q_LR(mpobw{k+1,j},mpobw{k+1,j+1},g2bw{j},sv_min,D);
    end
    %apply odd  2q terms
    for j=1:2:n-1
        [mpobw{k+1,j},mpobw{k+1,j+1}]=gate_2q_LR(mpobw{k+1,j},mpobw{k+1,j+1},g2bw{j},sv_min,D);
    end
    %apply 1q terms
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalisation
if mod(k-1,midstep)==0||k==nbin
mpofw(k+1,:)=mpo_normalize(mpofw(k+1,:));
mpobw(k+1,:)=mpo_normalize(mpobw(k+1,:));
end
maxD=1;
for j=1:n
    maxD=max(maxD,max(size(mpofw{k+1,j})));
end
maxDfw(k)=maxD;
maxD=1;
for j=1:n
    maxD=max(maxD,max(size(mpobw{k+1,j})));
end
maxDbw(k)=maxD;
normfw(k)=overlap(mpofw(k+1,:),mpofw(k+1,:));
normbw(k)=overlap(mpobw(k+1,:),mpobw(k+1,:));
end
for k=1:nbin
    ovl_check(k)=overlap(mpobw(nbin-k+2,:),mpofw(k,:));
end
end
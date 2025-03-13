function [iF,iG]=mps_infid(H0,Hc,c,T,mps0,mpstg,sv_min,D,midstep,nt)
n=length(mps0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
d=size(mps0{1},2);
dt=T/(nbin*nt);
Dt=T/nbin;
mpsfw=cell(nbin+1,n);  
mpsbw=cell(nbin+1,n);
iG=zeros(nbin,nc);
%%%%%%%%%%%%%%%%%%%%%%%%
mpsfw(1,:,:)=mps0;
mpsbw(1,:,:)=mpstg;
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
   mpsfw(k+1,:)=mpsfw(k,:);
    for jt=1:nt 
     %apply 1q terms
     for j=1:n
        [mpsfw{k+1,j}]=mps_gate_1q(mpsfw{k+1,j},g1{j});
    end
     %apply odd 2q terms
    for j=1:2:n-1
        [mpsfw{k+1,j},mpsfw{k+1,j+1}]=mps_gate_2q(mpsfw{k+1,j},mpsfw{k+1,j+1},g2{j},sv_min,D);
    end
     %apply even 2q terms
    for j=2:2:n-1
        [mpsfw{k+1,j},mpsfw{k+1,j+1}]=mps_gate_2q(mpsfw{k+1,j},mpsfw{k+1,j+1},g2{j},sv_min,D);
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
   mpsbw(k+1,:)=mpsbw(k,:);
    %left apply
    for jt=1:nt
    %apply even 2q terms
    for j=2:2:n-1
        [mpsbw{k+1,j},mpsbw{k+1,j+1}]=mps_gate_2q(mpsbw{k+1,j},mpsbw{k+1,j+1},g2bw{j},sv_min,D);
    end
    %apply odd  2q terms
    for j=1:2:n-1
        [mpsbw{k+1,j},mpsbw{k+1,j+1}]=mps_gate_2q(mpsbw{k+1,j},mpsbw{k+1,j+1},g2bw{j},sv_min,D);
    end
    %apply 1q terms
    for j=1:n
        [mpsbw{k+1,j}]=mps_gate_1q(mpsbw{k+1,j},g1{j});
    end
    end
    %normalisation
if mod(k-1,midstep)==0||k==nbin
mpsfw(k+1,:)=mps_normalize(mpsfw(k+1,:));
mpsbw(k+1,:)=mps_normalize(mpsbw(k+1,:));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gradient
ovl=mps_overlap(mpstg,mpsfw(nbin+1,:));
iF=1-abs(ovl)^2;
for k=1:nbin
    tnsfw_temp=mpsfw(k+1,:);
    tnsbw_temp=mpsbw(nbin-k+1,:);
%     ovl_check_temp=overlap_TN(tnsbw_temp,tnsfw_temp);
    for jc=1:nc  
            ovl_diff_left=0;
%             ovl_diff_right=0;
           for js=1:length(Hc(jc).sys)
                 jq=Hc(jc).sys(js);
                 gate=-1i*Dt*Hc(jc).op{js};
                 tnsfw_diff_left=tnsfw_temp;
%                  tnsfw_diff_right=conjtp(tnsfw_temp);
            tnsfw_diff_left{jq}=mps_gate_1q(tnsfw_diff_left{jq},gate);
%             tnsfw_diff_right{jm,jq}=gate_1q(tnsfw_diff_right{jm,jq},gate);
%              tnsfw_diff_right=conjtp(tnsfw_diff_right);
            ovl_diff_left=ovl_diff_left+mps_overlap(tnsbw_temp,tnsfw_diff_left);
%             ovl_diff_right=ovl_diff_right+overlap_TN(tnsbw_temp,tnsfw_diff_right);
           end
%              ovl_diff=ovl_diff_left+ovl_diff_right;
              iG(k,jc)=-2*real(ovl_diff_left*conj(ovl));
    end
end
iG=iG(:);
end

function [iF,iG]=infid_rcp(H0,Hc,c,T,mpo0,mpotg,sv_min,D,midstep,nt)
n=length(mpo0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
d=size(mpo0{1},2);
dt=T/(nbin*nt);
Dt=T/nbin;
mpofw=cell(nbin+1,n);  
mpobw=cell(nbin+1,n);
iG=zeros(nbin,nc);
%%%%%%%%%%%%%%%%%%%%%%%%
mpofw(1,:,:)=mpo0;
mpobw(1,:,:)=mpotg;
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
    %apply 1q terms
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
    end
    %normalisation
if mod(k-1,midstep)==0||k==nbin
mpofw(k+1,:)=mpo_recompress(mpofw(k+1,:),sv_min,D,2);
mpobw(k+1,:)=mpo_recompress(mpobw(k+1,:),sv_min,D,2);
mpofw(k+1,:)=mpo_normalize(mpofw(k+1,:));
mpobw(k+1,:)=mpo_normalize(mpobw(k+1,:));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gradient
for k=1:nbin
    tnsfw_temp=mpofw(k+1,:);
    tnsbw_temp=mpobw(nbin-k+1,:);
%     ovl_check_temp=overlap_TN(tnsbw_temp,tnsfw_temp);
    for jc=1:nc  
            ovl_diff_left=0;
%             ovl_diff_right=0;
           for js=1:length(Hc(jc).sys)
                 jq=Hc(jc).sys(js);
                 gate=-1i*Dt*Hc(jc).op{js};
                 tnsfw_diff_left=tnsfw_temp;
%                  tnsfw_diff_right=conjtp(tnsfw_temp);
            tnsfw_diff_left{jq}=gate_1q(tnsfw_diff_left{jq},gate);
%             tnsfw_diff_right{jm,jq}=gate_1q(tnsfw_diff_right{jm,jq},gate);
%              tnsfw_diff_right=conjtp(tnsfw_diff_right);
            ovl_diff_left=ovl_diff_left+overlap(tnsbw_temp,tnsfw_diff_left);
%             ovl_diff_right=ovl_diff_right+overlap_TN(tnsbw_temp,tnsfw_diff_right);
           end
%              ovl_diff=ovl_diff_left+ovl_diff_right;
              iG(k,jc)=-2*real(ovl_diff_left);
    end
end
ovl=overlap(mpotg,mpofw(nbin+1,:));
iF=1-real(ovl);
iG=iG(:);
end

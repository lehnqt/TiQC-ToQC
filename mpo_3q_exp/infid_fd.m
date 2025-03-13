function [iF,iG] = infid_fd(H0,Hc,c,T,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso)
iF=infid_nograd(H0,Hc,c,T,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
if nargout>1
dc=10^(-10);
iG=zeros(length(c),1);
parfor jc=1:length(c)
    c1=c;
    c1(jc)=c1(jc)+dc;
    iF1=infid_nograd(H0,Hc,c1,T,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
    iG(jc)=(iF1-iF)/dc;
end
end
end
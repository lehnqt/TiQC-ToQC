%fix Nan problem
function wf=expv_krylov(t,H,wi,numK)
if issparse(H)==0
    fprintf('warning: maxtrix is not sparse!\n');
end
m=min(numK,length(H));
wf=zeros(size(wi));
for i=1:size(wi,2)
w0=wi(:,i);
wL=zeros(length(w0),m);
alpha=zeros(1,m);
beta=zeros(1,m);
v1=w0;
v0=0;
beta(1)=0;
j=1;
wL(:,1)=v1;
%energy%
nanflag=0;
while j<m
w=H*v1;
alpha(j)=w'*v1;
w=w-alpha(j)*v1-beta(j)*v0;
beta(j+1)=norm(w);
v0=v1;
v1=w/beta(j+1);
wL(:,j+1)=v1;
if isnan(alpha(j))
    nanflag=1;
    break;
end
j=j+1;
end
w=H*v1;
alpha(m)=w'*v1;
if nanflag==1
    M=j-1;
else
    M=m;
end
T=diag(alpha(1:M),0)+diag(beta(2:M),1)+diag(beta(2:M),-1);
T=expm(-1i*t*T);
u=T(1:M,1);
wf(:,i)=wL(:,1:M)*u;
end
end

% N=11;
D=10;
precision=1e-5;
% A=1;
% T=2;
% nbin=100;
% dt=T/nbin;
% jflipped=5;
% magnetization in z-direction
% c=ones(nbin,N,2);
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2);
s_plus=(sx+1i*sy)/2;
s_minus=(sx-1i*sy)/2;
p0=(id+sz)/2;
p1=(id-sz)/2;
factor=@(j) (eq(j,1)+eq(j,N)+1)/2; 
oset0=cell(1,N);
for j=1:N
    oset0{1,j}=sz; 
end
oset1=cell(1,N);
for j=1:N
    oset1{1,j}=p0; 
end
oset2=cell(1,N);
for j=1:N
    oset2{1,j}=p1; 
end
oset3=cell(1,N);
for j=1:N
    oset3{1,j}=s_plus; 
end
oset4=cell(1,N);
for j=1:N
    oset4{1,j}=s_minus; 
end
% oset{1,jflipped}=sx;
% time evolution operator
I=reshape(id,[1,1,2,2]);
mpo_even=cell(1,N);
mpo_odd=cell(1,N);
for j=1:N
    mpo_even{j}=I; 
    mpo_odd{j}=I;
end
mps0=cell(1,N);
for j=1:N
% if j==jflipped 
%     state=[0; 1]; 
% else
    state=[1; 0];
% end
mps0{j}=reshape(state,[1,1,2]);
end
% time evolution
mpsGHZ=cell(1,N);

mpsGHZ{1}=zeros(1,2,2);
mpsGHZ{1}(1,1,1)=2^(-1/(2*N));
mpsGHZ{1}(1,2,2)=2^(-1/(2*N));

for j=2:N-1
mpsGHZ{j}=zeros(2,2,2);
mpsGHZ{j}(1,1,1)=2^(-1/(2*N));
mpsGHZ{j}(2,2,2)=2^(-1/(2*N));
end
mpsGHZ{N}=zeros(2,1,2);
mpsGHZ{N}(1,1,1)=2^(-1/(2*N));
mpsGHZ{N}(2,1,2)=2^(-1/(2*N));


mps=mps0;
ghzoverlap_MPS=zeros(nbin,1);
for k=1:nbin
    fprintf('step %d\n',k);
for j=1:(N-1) 
%     h=kron(s_plus,s_minus)+kron(s_minus,s_plus)+(c(k,j,1)*(kron(sx,id))+c(k,j,2)*(kron(sy,id)))*factor(j)+(c(k,j+1,1)*kron(id,sx)+c(k,j+1,2)*kron(id,sy))*factor(j+1);
%     h=kron(sz,sz)+A*((c(k,j,1)*(kron(sx,id))+c(k,j,2)*(kron(sy,id)))*factor(j)+(c(k,j+1,1)*kron(id,sx)+c(k,j+1,2)*kron(id,sy))*factor(j+1));
    h=J(j)*(kron(s_plus,s_minus)+kron(s_minus,s_plus))+A*(c(k,j,1)*kron(sx,id)+c(k,j,2)*kron(sy,id));
    w=expm(-1i*dt*h);
    w=reshape(w,[2,2,2,2]); w=permute(w,[1,3,2,4]); w=reshape(w,[4,4]);
    [U,S,V]=svd2(w);
    eta=size(S,1);
    U=U*sqrt(S); V=sqrt(S)*V;
    U=reshape(U,[2,2,eta]); U=permute(U,[4,3,2,1]);
    V=reshape(V,[eta,2,2]); V=permute(V,[1,4,3,2]);
    if mod(j,2)==1
    mpo_odd{j}=U;
    mpo_odd{j+1}=V;
    else
    mpo_even{j}=U;
    mpo_even{j+1}=V;
    end 
end
mps=reduceD(mps,mpo_even,D,precision);
mps=reduceD(mps,mpo_odd,D,precision);
% ghzoverlap_MPS(k)=(expectationvalue(mps,oset1)+expectationvalue(mps,oset2)+expectationvalue(mps,oset3)+expectationvalue(mps,oset4))/2;
% ghzoverlap_MPS(k)=expectationvalue(mps,oset0);
ghzoverlap_MPS(k)=overlap(mpsGHZ,mps);
end
% starting state (one spin flipped)
time_grid=dt*(1:nbin);
hold on;
plot(time_grid,real(ghzoverlap_MPS));
hold on;
plot(time_grid,imag(ghzoverlap_MPS));
hold on;
plot(time_grid,abs(ghzoverlap_MPS).^2);
hold off;

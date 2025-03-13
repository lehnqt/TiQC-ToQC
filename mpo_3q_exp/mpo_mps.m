function mps=mpo_mps(mpo,mps0)%apply mpo on mps: 1:left, 2:down, 3:up, 4:right
n=length(mpo);
mps=mps0;
%complex transpose of tns
A=mpo{1};
B=mps0{1};%this is conjugation, not complex transpose
tensors={A,B};
connects={[-1,-3,1,-4],[-2,1,-5]};
O=ncon(tensors,connects);
O=reshape(O,[size(A,1)*size(B,1),size(A,2),size(A,4),size(B,3)]);
for j=2:n
    A=mpo{j};
    B=mps0{j};
    tensors={O,A,B};
    connects={[-1,-2,1,2],[1,-3,3,-4],[2,3,-5]};
    O=ncon(tensors,connects);
    T=reshape(O,[size(O,1)*size(O,2),size(A,2)*size(A,4)*size(B,3)]);
    [Utemp,Stemp,Vtemp]=svd(T,'econ');
    Vhtemp=Vtemp';
    dtemp=length(Stemp);
    mps{j-1}=reshape(Utemp*sqrt(Stemp),[size(O,1),size(O,2),dtemp]);
    O=reshape(sqrt(Stemp)*Vhtemp,[dtemp,size(A,2),size(A,4),size(B,3)]);
end
mps{n}=reshape(O,[size(O,1),size(O,2),size(O,3)*size(O,4)]);
end
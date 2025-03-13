function ovl=overlap(mpo,mpo0)%calculate tr(U+U0), inddices order of tns and tns0: 1:left, 2:down, 3:up, 4:right
%complex transpose of tns
n=length(mpo);
mpo=conjtp(mpo);
ovl=0;
%complex transpose of tns
A=mpo0{1};
B=mpo{1};%this is conjugation, not complex transpose
tensors={A,B};
connects={[-1,1,2,-3],[-2,2,1,-4]};
O=ncon(tensors,connects);
for j=2:n
    A=mpo0{j};
    B=mpo{j};
    tensors={O,A,B};
    connects={[-1,-2,1,2],[1,3,4,-3],[2,4,3,-4]};
    O=ncon(tensors,connects);
end
ovl=ovl+O;
end
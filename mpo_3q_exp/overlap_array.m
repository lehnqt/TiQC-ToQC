function ovl=overlap_array(tns,tns0)%calculate tr(U+U0), inddices order of tns and tns0: 1:left, 2:down, 3:up, 4:right
%complex transpose of tns
[m,n]=size(tns);
[m0,~]=size(tns0);
tns=conjtp_array(tns);
ovl=0;
%complex transpose of tns
for k0=1:m0
    for k=1:m 
A=tns0{k0,1};
B=tns{k,1};%this is conjugation, not complex transpose
tensors={A,B};
connects={[-1,1,2,-3],[-2,2,1,-4]};
O=ncon(tensors,connects);
for j=2:n
    A=tns0{k0,j};
    B=tns{k,j};
    tensors={O,A,B};
    connects={[-1,-2,1,2],[1,3,4,-3],[2,4,3,-4]};
    O=ncon(tensors,connects);
end
ovl=ovl+O;
    end
end
end
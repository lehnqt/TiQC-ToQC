function O=mps_overlap(mps,mps0)
N=length(mps0);
A=mps0{1};
B=conj(mps{1});
tensors={A,B};
connects={[-1,1,-3],[-2,1,-4]};
O=ncon(tensors,connects);
for j=2:N
    A=mps0{j};
    B=conj(mps{j});
    tensors={O,A,B};
    connects={[-1,-2,1,2],[1,3,-3],[2,3,-4]};
    O=ncon(tensors,connects);
end
end
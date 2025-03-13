function H=hsb(J,EB,ispbc)
Sx = sparse([0, 1; 1, 0]/2);
Sy = sparse([0, -1i; 1i , 0]/2);
Sz = sparse([1, 0; 0, -1]/2);
Ns=length(J);
H=sparse(2^Ns,2^Ns);
for k=1:Ns-1
    H=H+J(k)*(kron(speye(2^(k-1)),kron(kron(Sx,Sx),speye(2^(Ns-k-1))))...
        +kron(speye(2^(k-1)),kron(kron(Sy,Sy),speye(2^(Ns-k-1))))...
       + kron(speye(2^(k-1)),kron(kron(Sz,Sz),speye(2^(Ns-k-1)))))...
       +EB(k)*kron(speye(2^(k-1)),kron(Sz,speye(2^(Ns-k))));
end
k=Ns;
H=H+EB(k)*kron(speye(2^(k-1)),kron(Sz,speye(2^(Ns-k))));
if ispbc==1
    H=H+J(k)*(kron(kron(Sx,speye(2^(Ns-2))),Sx)...
        + kron(kron(Sy,speye(2^(Ns-2))),Sy)...
       + kron(kron(Sz,speye(2^(Ns-2))),Sz));
end
end

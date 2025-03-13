%index order of A0: 1,2,3,4 == left, bottom,top, right
function A=gate_1q_LR(A0,U) %apply the operator U of a single qubit gate on single-qubit tensor A0
Uc=U';
tensors={A0,U,Uc};
connects={[-1,1,2,-4],[-2,1],[2,-3]};
A=ncon(tensors,connects);
end
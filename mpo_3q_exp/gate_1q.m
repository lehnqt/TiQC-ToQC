%index order of A0: 1,2,3,4 == left, bottom,top, right
function A=gate_1q(A0,U) %apply the operator U of a single qubit gate on single-qubit tensor A0
tensors={A0,U};
connects={[-1,1,-3,-4],[-2,1]};
A=ncon(tensors,connects);
end
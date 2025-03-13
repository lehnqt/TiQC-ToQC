function A=mps_gate_1q(A0,tno) %apply the tensor network operator tno of a single qubit gate on single-qubit tensor A0
tensors={A0,tno};
connects={[-1,1,-3],[-2,1]};
A=ncon(tensors,connects);
end
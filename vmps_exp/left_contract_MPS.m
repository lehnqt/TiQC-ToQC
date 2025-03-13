function [sigBA,sigAB] = left_contract_MPS(sigBA,sBA,A,sAB,B)
% Contract an infinite 2-site unit cell from the left for the environment
% density matrices sigBA (B-A link) and sigAB (A-B link)

% initialize the starting vector
chiBA = size(A,1);
if size(sigBA,1) == chiBA
  StartVector = sigBA(:);
else
  StartVector = reshape(eye(chiBA),[chiBA^2,1]);
end

% define network for transfer operator contraction
tensors = {diag(sBA),diag(sBA),A,conj(A),diag(sAB),diag(sAB),B,conj(B)};
labels = {[1,2],[1,3],[2,4],[3,5,6],[4,5,7],[6,8],[7,9],[8,10,-1],[9,10,-2]};

%define function for boundary contraction and pass to eigs
left_iter = @(sigBA,tensors,labels,chiBA)(reshape(ncon(...
  [reshape(sigBA,[chiBA,chiBA]),tensors],labels),[chiBA^2,1]));
[sigBA,~] = eigs(@(sigBA)left_iter(sigBA(:),tensors,labels,chiBA),...
  chiBA^2,1,'largestabs','Display',0,'IsFunctionSymmetric',0,...
  'MaxIterations',300,'StartVector',StartVector);

% normalize the environment density matrix sigBA
if isreal(A)
  sigBA = real(sigBA);
end
sigBA = reshape(sigBA,[chiBA,chiBA]);
sigBA = 0.5*(sigBA + sigBA')./trace(sigBA);

% compute density matric sigAB for A-B link
sigAB = ncon({sigBA,diag(sBA),diag(sBA),A,conj(A)},...
  {[1,2],[1,3],[2,4],[3,5,-1],[4,5,-2]});
sigAB = sigAB./trace(sigAB);
end
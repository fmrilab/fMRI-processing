function l = depl2(X,t_k,mask2d)
% function l = depl2(X,t_k,mask2d)
% This is the dependent log-likelihood model from the Genovese paper.
%   X = vector of unknowns to be determined
%   t_k = data array (counts)
%
% Jon-Fredrik Nielsen, jfnielse@umich.edu
% $Id: depl2.m,v 1.2 2013/11/05 12:53:51 jfnielse Exp $

K = size(t_k,3);
pA = [X(1:K)];
pI = [X((K+1):(2*K))];
lambda = X(2*K+1);

[I,J] = find(mask2d);
tk1 = ones(K,1);
l = 0;
for pi = 1:length(I)
	tk1 = t_k(I(pi),J(pi),:);
	l = l + log(lambda*prod(pA(:).^tk1(:)) + (1-lambda)*prod(pI(:).^tk1(:)));
end

l = -l;     

return;


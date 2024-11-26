function [c,ceq] = mycon(X)
% nonlinear inequality and equality constraints, use with fmincon in demo.m
%
% $Id: mycon.m,v 1.2 2013/11/05 12:53:51 jfnielse Exp $

K = (length(X)-1)/2;
ceq = [sum(X(1:K))-1 sum(X((K+1):(2*K)))-1];  % the probabilities must sum to 1
c = []; % c is the inequality constraint

return;

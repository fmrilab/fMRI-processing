function [pA,pI,lambda] = demo(CT,Xinit,cmaps,mask2d,lambda_upper)
% Calculate ROC curve using dependent likelihood model from Genovese et al, see:
% Estimating test-retest reliability in functional MR imaging. I: Statistical methodology.
% Genovese CR, Noll DC, Eddy WF.
% Magn Reson Med. 1997 Sep;38(3):497-507. Review.
% and
% Estimating test-retest reliability in functional MR imaging. II: Application to motor and cognitive activation studies.
% Noll DC, Genovese CR, Nystrom LE, Vazquez AL, Forman SD, Eddy WF, Cohen JD.
% Magn Reson Med. 1997 Sep;38(3):508-17.
%
% INPUT:
%   CT = array of activation map thresholds
%   Xinit = initial guess for the unknowns (see below)
%   cmap = statistical activation map (e.g., correlation or t-map) of size [NxNxMxK], where 
%      N = image size
%      M = number of fmri run repetitions
%      K = number of activation thresholds
% Out:
% pA: detection rate;
% pI: false positive rate;
% lambda: fraction of true activated pixel
%
% Jon-Fredrik Nielsen, jfnielse@umich.edu
% Modified by Hao Sun, sunhao@umich.edu
%
% $Id: demo.m,v 1.3 2013/11/05 12:53:51 jfnielse Exp $

if ~exist('CT','var')
	CT = [0.2:0.05:0.5]; % threshold list
   K = length(CT); 
   Xinit = [ones(1,K+1)/(K+1) ones(1,K+1)*0.02 0.05];
end
if ~exist('lambda_upper','var')
   lambda_upper = 0.1;
end

% load data and roi
if ~exist('cmaps','var')
   datdir = './data/'; 
	load(sprintf('%s/roi.mat',datdir));
	mask2d = roi;
	for ii = 1:length(CT)
		load(sprintf('%s/cmaps-%.2f-clustersize1.mat',datdir,CT(ii)));   % load 'cmap' structure
		cmaps(:,:,:,ii) = cmap.bold;
	end
end

M = size(cmaps,3);   % number of fmri replications
K = length(CT);      % number of thresholds

% apply spatial mask
for k = 1:K
	for m = 1:M
		cmap2d = cmaps(:,:,m,k);
		cmap2d(~mask2d) = NaN;
		cmaps(:,:,m,k) = cmap2d;
	end
end

% count
cmaps(abs(cmaps)>0) = 1;      % binarize activation maps to facilitate counting
cnt = squeeze(sum(cmaps,3));  % number of times (out of M) a voxel is classified as active, for each k (threshold level)
t_k(:,:,1) = M - cnt(:,:,1);  % t_0, number of times a voxel is not classified as active for any k
t_k(:,:,K+1) = cnt(:,:,K);    % t_K, number of times a voxel is classified active for all threshold levels
for k = 2:K
	t_k(:,:,k) = cnt(:,:,k-1)-cnt(:,:,k); % number of times a voxel is classified as active at (k-1)-th threshold level
end

% solve the MLE 
options = optimset('MaxFunEvals',1e5);
options = optimset(options,'Algorithm','sqp');
AA = zeros(size(Xinit)); AA(end) = 1; 
[X Fval] = fmincon(@(x) depl2(x,t_k,mask2d),Xinit',AA,lambda_upper,[],[],0*Xinit,ones(size(Xinit)),@mycon,options);

% calculate the propabilities pA(k) and pA(k), i.e., the estimated true/false positive rate for threshold k
pAd = [X(1:(K+1))];
pId = [X((K+2):(2*(K+1)))];
lambda = X(2*(K+1)+1);
for k = 2:(K+1)
	pA(k-1) = sum(pAd(k:(K+1)));
	pI(k-1) = sum(pId(k:(K+1)));
end

[pA pI lambda]
figure
plot(pI,pA,'bo-');
xlabel('false positive rate');
ylabel('detection rate');
title('ROC'); 
return;


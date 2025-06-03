function [tscore,beta] = fmri_tscore(A,X,df)
% function to calculate fMRI tscore maps from given timeseries and design
% matrix
% by David Frey
%
% inputs:
% A - design matrix (nt x nc)
% X - image timeseries ((N) x nt)
% df = degrees of freedom (optional, default = nt - nc)
%
% outputs:
% tscore - voxel-wise activation tscores ((N) x nc)
% beta - voxel-wise activation beta coefficients ((N) x nv)
%

    % get sizes
    nt = size(A,1); % number of time points
    nc = size(A,2); % number of contrasts
    xvec = reshape(X,[],nt).'; % vectorized timeseries
    sz = size(X);
    nv = size(xvec,2); % number of voxels

    % set default degrees of freedom
    if nargin<3 || isempty(df)
        df = nt - nc;
    end

    % calculate beta coefficients
    beta = pinv(A) * xvec;

    % calculate variance
    r = xvec - A*beta; % residual
    rss = sum(r.^2, 1); % sum of squares
    v = rss / df; % variance (sigma^2)

    % loop through voxels
    tscore = zeros(nc,nv);
    A_gram = A'*A;
    for n = 1:nv
        SE_beta = sqrt(v(n) * diag(pinv(A_gram))); % standard error of beta
        tscore(:,n) = beta(:,n) ./ SE_beta(:);
    end

    % reshape to image dimensions
    tscore = reshape(tscore.',[sz(1:end-1),nc]);
    beta = reshape(beta.',[sz(1:end-1),nc]);

end


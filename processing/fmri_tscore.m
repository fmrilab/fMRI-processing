function [tscore,beta] = fmri_tscore(A,y,df)
% A = [nt x nc] design matrix
% y = [nt x nv] vectorized image timeseries
% df = degrees of freedom (optional, default = nt - nc)
% tscore = [nc x nv] vectorized voxel-wise activation tscores
% beta = [nc x nv] vectorized voxel-wise activation beta coefficients

    % get sizes
    nt = size(A,1); % number of time points
    nc = size(A,2); % number of contrasts
    nv = size(y,2); % number of voxels

    % set default degrees of freedom
    if nargin<3 || isempty(df)
        df = nt - nc;
    end

    % calculate beta coefficients
    beta = pinv(A) * y;

    % calculate variance
    r = y - A*beta; % residual
    rss = sum(r.^2, 1); % sum of squares
    v = rss / df; % variance (sigma^2)

    % loop through voxels
    tscore = zeros(nc,nv);
    A_gram = A'*A;
    for n = 1:nv
        SE_beta = sqrt(v(n) * diag(pinv(A_gram))); % standard error of beta
        tscore(:,n) = beta(:,n) ./ SE_beta(:);
    end

end


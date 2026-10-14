function X_detrend = poly_detrend(X,n)
% function to remove polynomial trends from temporal data
% except for the DC component because presumably a DC regressor will be
% included in the design matrix in GLM
% by David Frey
%
% inputs:
% X - image timeseries ((N) x nt)
% n - max polynomial degree
%
% outputs:
% X_out - detrended data ((N) x nt)
%

    % create time vector
    t = (0:size(X,ndims(X))-1).';

    % create polynomial regressors
    A = t .^ (0:n);

    % get coefficients
    [~,coeffs] = fmri_tscore(A,X);
    coeffs = reshape(coeffs,[],n+1);

    % reshape into trends image
    trends = reshape(coeffs(:,2:end)*A(:,2:end).', size(X));

    % subtract out trends from magnitude
    X_detrend = X - trends;

end
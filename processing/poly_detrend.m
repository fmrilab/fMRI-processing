function X_detrend = poly_detrend(X,n)
% function to remove polynomial trends from temporal data
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
    [~,coeffs] = fmriutl.fmri_tscore(A,abs(X));
    coeffs = reshape(coeffs,[],n+1);

    % reshape into trends image
    trends = reshape(coeffs(:,2:end)*A(:,2:end).', size(X));

    % subtract out trends from magnitude
    X_detrend = (abs(X) - trends).*exp(1i*angle(X));

end
function s = GetRepeatability(X, nSample)

% GETREPEATABILITY compute a measure of state repeatability of evolution
%
% s = GetRepeatability(X, nSample) computes the euclidean-distance based
% repeatability index. It samples nSample replicates to evaluate the index
% at each time point
%
% INPUT
%
% X is a nT by nD by nRep array of relative abundnaces of nD species across
% nT time points in nRep replicates. nSample is the number of random
% replicates to sample for evaluating the repeatability index at any time
% point.
%
% OUTPUT
% 
% s is an nT by 1 vector of the repeatability index over time.

nT = size(X,1); % number of time points
nD = size(X,2); % number of dimensions of the measurements
nRep = size(X,3); % number of replicates

s = nan(nT,1);

for it=1:nT
    Y1 = squeeze(X(it,:,randsample(nRep, nSample)))';
    Y2 = squeeze(X(it,:,randsample(nRep, nSample)))';
    
    D = sqrt(sum((Y1 - Y2).^2,2));
    s(it) = 1 - mean( D ) / sqrt(2) ;    
end

% NNPC for precomputed dissimilarity matrix
% D: Dissimilarity matrix
% q: Number of nearest neighbors
% L: Number of clusters
% labels: Labels assigned to the data points


function labels = NNPCDist(D,q,L)
    N = size(D,1);

    Z = zeros(N,N);

    % Compute nearest neighbor graph
    for i=1:N
        di = D(:,i);
        di(i) = inf; % avoid that the ith point is selected
        [dio,order] = sort(di, 'ascend');
        Z(i, order(1:q) ) = exp(-2*dio(1:q));
    end

    Z = Z + Z';

    % Normalized spectral clustering
    Dsum = diag(1./sqrt(sum(Z)+eps));
    LN = eye(N) - Dsum*Z*Dsum;
    [~,S,V] = svd(LN);
    
    if nargin < 3
        singvals = diag(S);
        [~,minidx] = min(diff(singvals(1:(end-1))));
        L = N - minidx;
    end
    
    VL = V(:,N-L+1:N);
    VL = normr(VL);
    [labels,~] = kmeans(VL,L,'maxiter',1000,'replicates',200,'EmptyAction','singleton');
end



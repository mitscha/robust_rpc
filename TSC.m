% TSC algorithm
% X: Data matrix containing data points as columns
% q: Number of nearest neighbors
% L: Number of clusters
% labels: Labels assigned to the data points


function labels = TSC(X,q,L)
    N = size(X,2);

    Z = zeros(N,N);
    
    D = abs(X'*X);
    Xnorm = normc(X);
    Dnorm = abs(Xnorm'*Xnorm);

    % Compute nearest neighbor graph
    for i=1:N
        di = D(:,i);
        di(i) = -inf; % avoid that the ith point is selected
        [~,order] = sort(di, 'descend');
        dion = Dnorm(order,i);
        Z(i, order(1:q) ) = exp(-2*acos(dion(1:q)));
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



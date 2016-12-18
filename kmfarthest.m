% k-means with furthest point initialization
% D: Dissimilarity matrix
% L: Number of clusters
% labels: Labels assigned to the data points

function labels = kmfarthest(D,L,varargin)
    idx = zeros(1,L);
    idx(1) = 1;
    
    for i = 1:L-1
        [~,idx(i+1)] = max(min(D(:,idx(1:i)),[],2));
    end
    
    if ~isempty(varargin)
        X = varargin{1};
        labels = kmeans(X',L,'Distance','cityblock','Start',X(:,idx)');
    else
        [~,labels] = min(D(:,idx),[],2);
    end
end
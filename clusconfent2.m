% Compute the entropy of the 2x2 clustering confustion matrix corresponding
% to labels 'labels' and the ground truth 'gt' (column or row vectors)

function S = clusconfent2(labels,gt)
    % Swap label indices if needed
    if sum(labels ~= gt) > 0.5
        labels = labels - 1;
        labels(labels == 0) = 2;
    end
    
    % Compute clustering confusion matrix
    M = zeros(2,2);
    
    M(1,1) = sum(labels(gt == 1) == 1);
    M(1,2) = sum(labels(gt == 2) == 1);
    M(2,1) = sum(labels(gt == 1) == 2);
    M(2,2) = sum(labels(gt == 2) == 2);
    
    % Compute entropy
    S = 0;
    for i = 1:2
        for j = 1:2
            if M(i,j) ~= 0
                S = S + M(i,j)*log(sum(M(i,:))/M(i,j));
            end
        end
    end
    S = S/sum(M(:));
end
% Compute minimum clustering error w.r.t. all label permutations (brute
% force)
% labels: Vector containing label assignments (values in 1:L)
% gt: Vector containing ground truth labels (values in 1:L)

function cemin = computece(labels,gt)
    % Transpose if needed
    if size(labels,2) < size(labels,1); labels = labels'; end
    if size(gt,2) > size(gt,1); gt = gt'; end
    
    L = max(gt);
    N = length(gt);
    cemin = 1;
    % L! label permutations
    per = perms(1:L)';
    
    indmatrix = (per(:,1)*ones(1,N) == ones(L,1)*labels);
    for i = 1:size(per,2)
        % Map labels and compute clustering error
        mapmatrix = repmat(per(:,i),1,N);
        ce = sum(mapmatrix(indmatrix) ~= gt)/N;
        if ce < cemin
            cemin = ce;
        end
    end
end
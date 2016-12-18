% Function to reproduce Figure 4 and Table 2 in "Robust nonparametric 
% nearest neighbor random process clustering" by Michael Tschannen and 
% Helmut Boelcskei

function testeeg

    s = RandStream('mcg16807','Seed',1);
    RandStream.setGlobalStream(s);

    infolders{1} = './S/'; % Set E: Seizure activity
    infolders{2} = './Z/'; % Set A: Healthy volunteers
    
    qs = 2:1:20;
    L = 2;

    winwidths = (20:20:1000)-1;
    
    X = [];
    gt = [];

    % Load data
    for i = 1:length(infolders)
        curfolder = dir(infolders{i});
        for j = 3:length(curfolder)
            xc = dlmread(strcat(infolders{i},curfolder(j).name));
            X = [X xc];
            gt = [gt i];
        end
    end
    
    % Subtract mean
    X = X - ones(size(X,1),1)*mean(X,1);
    
    
    
    % Shuffle observations
    idxperm = randperm(size(X,2));
    X = X(:,idxperm);
    gt = gt(idxperm);
    
    len = size(X,1);
    
    cesnnpc = zeros(length(winwidths),length(qs));
    tcompsnnpc = cesnnpc;
    ceskm = zeros(length(winwidths),1);
    tcompskm = ceskm;
    ceskmfp = ceskm;
    tcompskmfp = ceskm;
    
    % Evaluate algorithms
    for i = 1:length(winwidths)
        w = zeros(2*len-1,1);
        w((len-(winwidths(i)-1)/2):(len+(winwidths(i)-1)/2)) = bartlett(winwidths(i));

        for j = 1:length(qs)
            if j == 1
                cesnnpc(i,j) = procclus(X,gt,w,qs(j),L,false);
            else
                [ce, tcomp] = procclus(X,gt,w,qs(j),L,true);
                cesnnpc(i,j) = ce(1);
                ceskm(i) = ce(2);
                ceskmfp(i) = ce(3);
                
                tcompsnnpc(i,j) = tcomp(1);
                tcompskm(i) = tcomp(2);
                tcompskmfp(i) = tcomp(3);
            end
        end
    end
    
    [mincennpc, minidxnnpc] = min(cesnnpc(:));
    [idx,idy] = ind2sub(size(cesnnpc),minidxnnpc);
    minq = qs(idy);
    
    [mincekm,minidxkm] = min(ceskm);
    [mincekmfp,minidxkmfp] = min(ceskmfp);
    
    fprintf('NNPC: CE = %f, q = %i, winwidth = %i, t = %f \n', mincennpc, minq, winwidths(idx),tcompsnnpc(idx,idy))
    fprintf('KM: CE = %f, winwidth = %i, t = %f \n', mincekm, winwidths(minidxkm),tcompskm(minidxkm))
    fprintf('KMfp: CE = %f, winwidth = %i, t = %f \n', mincekmfp, winwidths(minidxkmfp),tcompskmfp(minidxkmfp))
    
    cesnnpc = permute(squeeze(cesnnpc),[2 1 3]);
    
    saveheatmap(cesnnpc,qs,winwidths,'cesnnpc.dat')
    
    dlmwrite('ceskm.dat',[winwidths' ceskm ceskmfp],'delimiter',' ','precision',4)
    
end

% Function evaluating the three algorithms for a given set of parameters
function [ce, tcomp] = procclus(X,gt,w,q,L,km)
    
    tic
    Xf = zeros(2*size(X,1),size(X,2));
    for k = 1:size(X,2)
        Xf(:,k) = abs(fft([0; w.*xcorr(X(:,k))]));
        % abs is needed here as xcorr seems to introduce a minor phase in
        % the time domain, which then results in high frequency modulation
        % in the frequency domain
    end

    Xf = Xf*diag(1./sum(Xf,1));

    D = zeros(size(Xf,2));
    for k = 1:size(Xf,2)
        for l = k:size(Xf,2)
            D(k,l) = sum(abs(Xf(:,k)-Xf(:,l)))/2;
            D(l,k) = D(k,l);
        end
    end
    tDist = toc;

    tic
    labelsNNPC = NNPCDist(D,q,L);
    tcomp = tDist + toc;
    
    ce = computece(labelsNNPC,gt);
    
    if km
        tic
        labelsKM = kmfarthest(D,L);
        tKM = tDist + toc;
        tic
        labelsKMfp = kmfarthest(D,L,Xf);
        tKMfp = tDist + toc;
        
        ce = [ce computece(labelsKM,gt) computece(labelsKMfp,gt)];
        tcomp = [tcomp tKM tKMfp];
    end
end
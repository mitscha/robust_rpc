% Function to reproduce Table 1 in "Robust nonparametric nearest neighbor 
% random process clustering" by Michael Tschannen and Helmut Boelcskei

% mode: 1: Subject #16, 2: Subject #35

function testMOCAP(mode)
    s = RandStream('mcg16807','Seed',1);
    RandStream.setGlobalStream(s);
    
    % Index of the marker sequence used for clustering
    markeridx = 53;
    
    L = 2;
    q = 5;

    if mode == 2
        infolder = '../MOCAP_S35/';
        N = 33;
        % Maximum observation length
        maxlen = 497;
        % Ground truth
        gt = ones(N,1);
        gt(17:26) = 2;
        
        exidxw = -1;
        exidxr = -1;
    else
        infolder = './MOCAP_S16/';
        N = 49;
        maxlen = 580;
        gt = ones(N,1);
        gt([1 26:37 39:48]) = 2;
        
        exidxw = 47;
        exidxr = 7;
    end
    
    files = dir(infolder);

    X = zeros(2*maxlen-1,N);
    
    Y = zeros(maxlen,N);

    % Offset file list index to skip '.', '..' and hidden files
    startidx = length(files) - N + 1;
    

    tDist = 0;
    for i = startidx:length(files)
        data = amc_to_matrix(strcat(infolder, files(i).name));
        len = length(data(:,1));
        
        tic
        w = bartlett(2*len-1);
        
        % Estimate PSD via Blackman-Tukey PSD estimator
        x = abs(fft([0; w.*xcorr(data(:,markeridx))],2*maxlen-1));
        X(:,i-startidx+1) = x;
        tDist = tDist + toc;
        
        Y(1:len,i-startidx+1) = data(:,markeridx);
    end

    % Shuffle data
    idxperm = randperm(N);
    X = X(:,idxperm);
    gt = gt(idxperm);

    tic
    % Normalixe PSD estimates to unit power
    X = X*diag(1./sum(X,1));

    % Compute dissimilarity matrix
    D = zeros(size(X,2));
    for k = 1:size(X,2)
        for l = k:size(X,2)
            D(k,l) = sum(abs(X(:,k)-X(:,l)))/2;
            D(l,k) = D(k,l);
        end
    end
    tDist = tDist + toc;

    tic
    labelsNNPC = NNPCDist(D,q,L);
    tNNPC = tDist + toc;
    tic
    labelsKM = kmfarthest(D,L);
    tKM = tDist + toc;
    tic
    labelsKMfp = kmfarthest(D,L,X);
    tKMfp = tDist + toc;
    fprintf('NNPC: CE = %f, S = %f, t = %f \nKM:   CE = %f, S = %f, t = %f \nKMfp: CE = %f, S = %f, t = %f \n ',...
        computece(labelsNNPC,gt),clusconfent2(labelsNNPC,gt),tNNPC,...
        computece(labelsKM,gt),clusconfent2(labelsKM,gt),tKM,...
        computece(labelsKMfp,gt),clusconfent2(labelsKMfp,gt),tKMfp)
end



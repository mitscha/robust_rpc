% Function to compare the clustering performance of NNPC and TSC, see
% Sec. 4 in "Robust nonparametric nearest neighbor random process
% clustering" by Michael Tschannen and Helmut Boelcskei

% sigmansq: Noise variance

function synthDataNNPCvsTSC(sigmansq)
    s = RandStream('mcg16807','Seed',1);
    RandStream.setGlobalStream(s);

    L = 3;
    nl = 25;
    q = 10;

    % Length of the realization
    m = 10000;
    % Observation lengths
    nlengths = 10;
    obslens = 2*floor(linspace(50, 2500, nlengths)/2);
    % Number of problem instances
    nrep = 20;


    % Width (variance) of the Gaussian window
    sigmawin = 50;

    % ARMA model parameter vectors
    a1 = 1;
    b1 = [0.75 1.0 -1.75 0.5];
    b1 = b1/sqrt(sum(b1.^2));

    a2 = 1;
    b2 = [0.5 1.25 -1.5 0.75];
    b2 = b2/sqrt(sum(b2.^2));

    a3 = [1 -0.2 0.4 0.1];
    a3 = a3*sqrt(sum(a3.^2));
    b3 = 1;
    
    % Plot model PSDs
    fab = @(a,b,n) abs(fft(b,n)./fft(a,n)).^2;

    nint = 1000;
    f1 = fab(a1,b1,nint);
    f2 = fab(a2,b2,nint);
    f3 = fab(a3,b3,nint);
    
    figure(1)
    plot([f1' f2' f3'])
    legend('s1(f)','s2(f)','s3(f)')
    xlabel('f')

    % Define ground truth
    gt = [ones(1,nl) 2*ones(1,nl) 3*ones(1,nl)];

    err = zeros(2,nlengths);

    for i = 1:nlengths
        curravg = zeros(1,2);
        for j = 1:nrep

            % Sample realizations
            X1 = filter(b1,a1,randn(m,nl));
            X2 = filter(b2,a2,randn(m,nl));
            X3 = filter(b3,a3,randn(m,nl));

            % Add noise
            if sigmansq > 0
                X1 = X1 + sqrt(sigmansq)*randn(m,nl);
                X2 = X2 + sqrt(sigmansq)*randn(m,nl);
                X3 = X3 + sqrt(sigmansq)*randn(m,nl);
            end

            XD = [X1 X2 X3];

            % Shuffle observations
            idxperm = randperm(size(XD,2));
            XD = XD(:,idxperm);
            gtcurr = gt(idxperm);

            % Estimate PSDs via Blackman-Tukey PSD estimator
            X = zeros(2*obslens(i),size(XD,2));
            for k = 1:size(XD,2)
                X(:,k) = abs(fft([0; gausswin(2*obslens(i)-1,(2*obslens(i)-1)/(2*sigmawin)).*xcorr(XD(end-obslens(i)+1:end,k))]));
            end

            % Normalize PSD estimates to unit power
            X = X*diag(1./sum(X,1));

            % Compute NNPC dissimilarity matrix
            D = zeros(size(X,2));
            for k = 1:size(X,2)
                for l = k:size(X,2)
                    D(k,l) = sum(abs(X(:,k)-X(:,l)))/2;
                    D(l,k) = D(k,l);
                end
            end
            

            labelsTSC = TSC(XD,q,L);
            curravg(1) = curravg(1) + computece(labelsTSC,gtcurr);

            labelsNNPC = NNPCDist(D,q,L);
            curravg(2) = curravg(2) + computece(labelsNNPC,gtcurr);

        end
        err(:,i) = curravg'/nrep;
    end
    
    % Plot CEs
    figure(2)
    plot(err')
    legend('TSC','NNPC')
    xlabel('M')
    ylabel('CE')
end

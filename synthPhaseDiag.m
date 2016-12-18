% Function to reproduce Figures 2 and 3 in "Robust nonparametric nearest 
% neighbor random process clustering" by Michael Tschannen and Helmut
% Boelcskei

function synthPhaseDiag

    s = RandStream('mcg16807','Seed',1);
    RandStream.setGlobalStream(s);

    nrep = 10;
    nl = 25;
    q = 10;
    m = 10000;
    winwidth = 101;
    
    nint = 1000;
    
    % Example PSDs plot
    s1 = fab(aresonator(0.4*pi, 0.6),1,nint);
    s2 = fab(aresonator(0.4*pi, 0.8),1,nint);
    s3 = fab(aresonator(0.7*pi, 0.6),1,nint);
    s4 = fab(aresonator(0.7*pi, 0.8),1,nint);
    
    dlmwrite(strcat('psds.dat'),[(0:(nint-1))'/nint s1'/sum(s1/nint) s2'/sum(s2/nint) s3'/sum(s3/nint) s4'/sum(s4/nint)],'delimiter',' ','precision',4)

    r = 0.6;
    phi0 = 0.7*pi;
    
    nint = 5000;
   
    % CE(obslength,distance) 
    suffix = '-distobslen';
    nlengths = 20;
    obslen = 2*floor(linspace(50, 2000, nlengths)/2);
    sigman = 0.5;
    phi = pi*(0.3:0.015:0.7);

    [ces,dX1X2] = testalgs(m,nl,q,phi0,r,phi,obslen,sigman,1,nrep,nint,winwidth);
    
    ces = permute(squeeze(ces),[2 1 3]);
    
    saveces(ces,obslen,dX1X2,suffix)   
    
    % CE(sigman,distance)
    suffix = '-distsigman';
    obslen = 300;
    sigman = 0:0.25:4;
    phi = pi*(0.1:0.015:0.7);
    
    [ces,dX1X2] = testalgs(m,nl,q,phi0,r,phi,obslen,sigman,1,nrep,nint,winwidth);

    ces = permute(squeeze(ces),[2 1 3]);
    
    saveces(ces,sigman,dX1X2,suffix) 
    
    % CE(obslength,distance)
    suffix = '-obslensigman';
    nlengths = 20;
    obslen = 2*floor(linspace(50, 2000, nlengths)/2);
    sigman = 0:0.25:4;
    phi = 0.62*pi;
    
    [ces,~] = testalgs(m,nl,q,phi0,r,phi,obslen,sigman,1,nrep,nint,winwidth);
    
    ces = squeeze(ces);
    
    saveces(ces,obslen,sigman,suffix)
    
    % CE(obslength,samplingprobability)
    suffix = '-obslenp';
    nlengths = 20;
    obslen = 2*floor(linspace(50, 2000, nlengths)/2);
    sigman = 0.5;
    phi = 0.62*pi;
    pminus = 1:0.3:7.2;
    
    [ces,~] = testalgs(m,nl,q,phi0,r,phi,obslen,sigman,1./pminus,nrep,nint,winwidth);
    
    ces = squeeze(ces);
    
    saveces(ces,obslen,pminus,suffix)    
end


function y = fab(a,b,n) 
    y = abs(fft(b,n)./fft(a,n)).^2;
end

function y = aresonator(phi,r) 
    y = [1 -2*r*cos(phi) r^2];
end


function [ces,dX1X2] = testalgs(m,nl,q,phi0,r,phi,obslen,...
    sigman,p,nrep,nint,winwidth)
    
    a1 = aresonator(phi0,r);
    b1 = 1;
    f1 = fab(a1,b1,nint);
    f1sum = sum(f1)/nint;
    b1 = 1/sqrt(f1sum);
    f1 = f1/f1sum;
    
    L = 2;
    
    ces = zeros(length(phi),length(obslen),length(sigman),length(p),3);
    dX1X2 = zeros(1,length(phi));
    
    % Define ground truth
    gt = [ones(1,nl) 2*ones(1,nl)];
    
    for i = 1:length(phi)
        
        a2 = [1 -2*r*cos(phi(i)) r^2];
        b2 = 1;
        f2 = fab(a2,b2,nint);
        f2sum = sum(f2)/nint;
        b2 = 1/sqrt(f2sum);
        f2 = f2/f2sum;
        
        dX1X2(i) = sum(abs(f1-f2)/2)/nint;
        for ip = 1:length(p)
            for j = 1:length(obslen)

                if obslen(j) >= winwidth
                    w = zeros(2*obslen(j)-1,1);
                    w((obslen(j)-(winwidth-1)/2):(obslen(j)+(winwidth-1)/2)) = bartlett(winwidth)/p(ip).^2;
                else
                    w = bartlett(winwidth)/p(ip).^2;
                    w = w(((winwidth+1)/2-obslen(j)+1):((winwidth+1)/2+obslen(j)-1));
                end
                w(obslen(j)) = w(obslen(j))*p(ip);

                for s = 1:length(sigman)

                    curravg = zeros(1,3);
                    for n = 1:nrep
                        % Sample realizations
                        X1 = filter(b1,a1,randn(m,nl));
                        X2 = filter(b2,a2,randn(m,nl));

                        % Add noise
                        if sigman(s) > 0
                            X1 = X1 + sigman(s)*randn(m,nl);
                            X2 = X2 + sigman(s)*randn(m,nl);
                        end

                        % Undersampling
                        if p(ip) < 1
                            X1 = (rand(m,nl) <= p(ip)).*X1;
                            X2 = (rand(m,nl) <= p(ip)).*X2;
                        end

                        XD = [X1 X2];

                        % Shuffle observations
                        idxperm = randperm(size(XD,2));
                        XD = XD(:,idxperm);
                        gtcurr = gt(idxperm);

                        % Estimate PSDs via Blackman-Tukey PSD estimator
                        X = zeros(2*obslen(j),size(XD,2));
                        for k = 1:size(XD,2)
                            % Need normalization here as xcorr computes the
                            % unnormalized autocorrelation
                            X(:,k) = abs(fft([0; w.*xcorr(XD(end-obslen(j)+1:end,k))]))/obslen(j);
                        end

                        % Compute dissimilarity matrix
                        D = zeros(size(X,2));
                        for k = 1:size(X,2)
                            for l = k:size(X,2)
                                D(k,l) = sum(abs(X(:,k)-X(:,l)))/(4*obslen(j));
                                D(l,k) = D(k,l);
                            end
                        end

                        labelsKM = kmfarthest(D,L);
                        curravg(1) = curravg(1) + computece(labelsKM,gtcurr);

                        labelsKMcont = kmfarthest(D,L,X);
                        curravg(2) = curravg(2) + computece(labelsKMcont,gtcurr);

                        labelsNNPC = NNPCDist(D,q,L);
                        curravg(3) = curravg(3) + computece(labelsNNPC,gtcurr);
                    end
                    ces(i,j,s,ip,:) = curravg/nrep;
                end
            end
        end
    end
end

% Save clustering error heatmaps
function saveces(ces,labx,laby,suffix)    
    saveheatmap(squeeze(ces(:,:,3)),labx,laby,strcat('cesnnpc',suffix,'.dat'));
    saveheatmap(squeeze(ces(:,:,1)),labx,laby,strcat('ceskm',suffix,'.dat'));
    saveheatmap(squeeze(ces(:,:,2)),labx,laby,strcat('ceskmfp',suffix,'.dat'));
end
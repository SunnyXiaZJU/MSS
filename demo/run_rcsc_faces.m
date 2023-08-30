
clear; clc;

load ../data/YaleB/YaleBCrop025.mat

lambda = 5;
alpha0 = 20;
beta = 5;
type = '1';

omega_iter = 10;
tau = 0.5;
epsilon = 1e-5;

opts.outiter = 100; 
opts.initer = 5; 
opts.gamma=1.1;
opts.orthflag = 0;
d = 9;

KSet = [2 3 5 8 10];
missrateTot = cell(max(KSet),1);
avgmissrate = -ones(max(KSet),1);
medmissrate = -ones(max(KSet),1);
tsTot = cell(numel(KSet),1);

missrateTot0 = cell(max(KSet),1);
avgmissrate0 = -ones(max(KSet),1);
medmissrate0 = -ones(max(KSet),1);
tsTot0 = cell(numel(KSet),1);


%parpool(4);
for i = 1:length(KSet)
    K = KSet(i);
    idx = Ind{K};  
    
    n = K*64;
    
    all_missrates = -ones(size(idx,1),1);
    ts = all_missrates; 
    
    all_missrates0 = -ones(size(idx,1),1);
    ts0 = all_missrates; 
    
    grps0 = s{K}';
    
    for j = 1:size(idx,1)    
        X = [];
        for p = 1:K
            X = [X Y(:,:,idx(j,p))];
        end
      
        opts.C = eye(n);
        Omega = ones(n)-eye(n);
        X = X*diag(1./sqrt(sum(X.*X)));
        alpha = alpha0/norm(X,1);
        tic
        
        for iter = 1 : omega_iter
            Omega_old = Omega;
            [C, out] = PRRSS(X, Omega, d*K, lambda, alpha, beta, type, opts);%rcsc_lnorm2m
            G = post_proc(C, 1, 0,1);
            [Omega, grps] = graph2omega(G,K,tau);
            
            if norm(Omega(:) - Omega_old(:), 1) < epsilon; break; end

            if iter == 1
                ts0(j) = toc;
                err = Misclassification(grps, grps0);
                all_missrates0(j) = err;
            end
        end
        ts(j) = toc;
        err = Misclassification(grps, grps0);
        
        all_missrates(j) = err;
        fprintf('data %d : missrate = %.4f, iter = %d \n', ...
            j, err, size(out.fvals,1));
        
       
    end
    tsTot{K} = ts;
    missrateTot{K} = all_missrates;
    avgmissrate(K) = mean(missrateTot{K});
    medmissrate(K) = median(missrateTot{K});
  
    disp([avgmissrate, medmissrate]);
    
    tsTot0{K} = ts0;
    missrateTot0{K} = all_missrates0;
    avgmissrate0(K) = mean(missrateTot0{K});
    medmissrate0(K) = median(missrateTot0{K});
end
save('out/faces/MSS_RO.mat', 'missrateTot', 'avgmissrate', 'medmissrate', 'tsTot', ...
    'missrateTot0', 'avgmissrate0', 'medmissrate0', 'tsTot0', 'opts');


% this code implements the algorithm for IEEE signal processing letter paper
%'Robust Subspace Clustering via Smoothed Rank Approximation'
% By Zhao Kang, 08/2015, Zhao.Kang@siu.edu
% This is for motion segmentation experiment with HopKins 155 data
clear all
close all
data = load_motion_data(1);

rho=67;%  Two tunable parameters
betazero=.4;

type=2;% different modeling of E
alpha=2;  % W cconstrct
iter=100;

errs = zeros(length(data), 1);
costs = zeros(length(data), 1);

for i = 1 : length(data)
    X = data(i).X;
    gnd = data(i).ids;
    K = max(gnd);
    [m,n]=size(X);
    
    %initializations
    W=eye(n);
    ww=ones(n,1);
    Y1=zeros(m,n);
    Y2=zeros(n);
    E=zeros(m,n);
    J=W;
    
    if abs(K - 2) > 0.1 && abs(K - 3) > 0.1
        id = i; % the discarded sequqnce
    end
    beta=betazero;
    for ii=1:iter
        if ii>=2
            Zold=Z;
        end
        
        
        gamma=W-Y2/beta;
        [ J] = LogSquare(gamma,beta/2);
        
        [E]=errormin(Y1,X,W,rho,beta,type);
        W=inv(eye(size(J))+X'*X)*(X'*(X-E)+J+(X'*Y1+Y2)/beta);
        Y1=Y1+beta*(X-X*W-E);
        Y2=Y2+beta*(J-W);
        beta=beta*1.1;
        Z=W;
        if type==1
            enorm=sum(sum(abs(X-X*Z)));
        elseif type==2
            enorm=sum(sum((X-X*Z).^2));
        else
            enorm=sum(sqrt(sum((X-X*Z).^2,1)));
        end
        func(ii)=log(det(eye(n)+Z'*Z))+rho*enorm;
        if ii>3 && norm(Z-Zold,'fro')/norm(Zold,'fro')<1e-5
            break
        end
        
    end
    
    if 1
    [U s V] = svd(Z);
    s = diag(s);
    r = sum(s>1e-6);
    
    U = U(:, 1 : r);
    s = diag(s(1 : r));
    V = V(:, 1 : r);
    
    M = U * s.^(1/2);
    mm = normr(M);
    rs = mm * mm';
    L = rs.^(2 * alpha);
    
    
    actual_ids = spectral_clustering(L, K);
    
    err = calculate_err(gnd, actual_ids);
    
    else
        grps = SpectralClustering(abs(Z)+abs(Z'), K);
        err = Misclassification(grps, gnd');
    end
    
    
    errs(i) = err;
    
    LK(i)=K;
    err2obj = max(0,mean(errs(LK==2)));
    mederr2obj = max(0,median(errs(LK==2)));
    err3obj = max(0,mean(errs(LK==3)));
    mederr3obj = max(0,median(errs(LK==3)));
    
    
    disp(['i:' num2str(i) ,',seg2 err=' num2str(err2obj) ',seg3 err=' num2str(err3obj) ',seg err=' num2str(err) ',rho=' num2str(rho) ', beta=' num2str(beta)]);
end

save('out/hpk/logdet.mat','errs','LK');
disp('results of all 156 sequences:');

disp('results of all motions:');
disp(['max = ' num2str(max(errs)) ',min=' num2str(min(errs)) ...
    ',median=' num2str(median(errs)) ',mean=' num2str(mean(errs)) ...
    ',std=' num2str(std(errs))] );

dlmwrite('all156result.txt', [rho betazero alpha ...
    max(errs) min(errs) median(errs) mean(errs) std(errs)], ...
    '-append', 'delimiter', '\t', 'newline', 'pc');


errs = errs([1:id-1,id+1:end]);
costs = costs([1:id-1,id+1:end]);
disp('results of all 155 sequences:');
disp(['max = ' num2str(max(errs)) ',min=' num2str(min(errs)) ...
    ',median=' num2str(median(errs)) ',mean=' num2str(mean(errs)) ',std=' num2str(std(errs))] );
dlmwrite('all155result.txt', [rho betazero alpha err2obj mederr2obj err3obj mederr3obj median(errs) mean(errs) std(errs) max(errs) min(errs)] , '-append', 'delimiter', '\t', 'newline', 'pc');



% final code for motion segmention

lambda = 10;
alpha = 50; 
beta = 0.05; 


tau = 0.5;
iniflag = 0;
epsilon = 0.1;
omega_iter = 10;


opts.gamma=1.1;
opts.epsilon1 = 1e-5;
opts.epsilon2 = 1e-5;
opts.outiter = 100;
opts.initer = 5;
opts.orthflag = 0; 
if exist('orthflag', 'var'),  opts.orthflag = orthflag; end

mu=2;
s0 = 1e-3;

d = 4;


data = load_motion_data(1);

errs = zeros(length(data), 1);
LK = zeros(length(data),1);
ts=  zeros(length(data),1);
errs0=  zeros(length(data),1);
id=-1;

for i = 1: length(data)
    fprintf('==== data = %d ===\n', i)
    X = data(i).X;
    gnd = data(i).ids;
    K = max(gnd);
    [~, n]=size(X);
    
    if K == 5 %abs(K - 2) > 0.1 && abs(K - 3) > 0.1
        id = i; % the discarded sequqnce
    end
    
    opts.C = eye(n);
    Omega = ones(n) - eye(n);
    
    tic;
    for iter = 1 : omega_iter
        Omega_old = Omega;
        [C, out] = PRRSS(X, Omega, d*K, lambda*n/K, alpha*n/K, beta*n/K, 'fro', opts);
        Z = BuildAdjacency(thrC(C,0.8));
        [U, s, ~] = svd(Z);
        s = diag(s); r = sum(s>s0);
        U = U(:, 1 : r); s = diag(s(1 : r));
        M = U * s.^(1/2); mm = normr(M);
        rs = mm * mm'; L = rs.^(2 * mu);
        [Omega, grps] = graph2omega(L, K, tau);
        if iniflag; opts.C = C; end
        temp1 = norm(Omega(:) - Omega_old(:), 1);
        if temp1 < epsilon; break; end
        err = Misclassification(grps, gnd);
        if iter == 1; errs0(i)=err;end
        fprintf('%d %.4f \n', length(out.fvals), err);        
    end
    ts(i) = toc;
    err = Misclassification(gnd, grps);
    errs(i) = err;
    
     fprintf('%d %.4f \n', length(out.fvals), err);
    
    LK(i)=K;
   
   
end

save('out/hpk/MSS_RO.mat', 'errs','errs0', 'LK', 'ts','id');

disp('results of all 155 sequences:');
motion_disp(errs, LK, id);
motion_disp(errs0, LK, id);

function [X, labels, settings] = gen_data(dk, nk, m, K, options, fname)
% generalize synthetic data for experiment in Section 8.2
% options : parameter besides dk, nk, K, namely r, randomflag, orthflag
% dk, nk : column vector or scalar represent dimension and size of each cluster
% K : the number of subspaces.
% m : length of one sample.
% X : each column is a sample.
% labels :  row column, each element from 1 to K.
% settings : parameters required to describe the data, which is different
% this script is the "general_inter" part of gen_data_All.m

% 10/04/2019 by Xia Yuqing

% set parameters
if numel(dk) == 1; dk = dk*ones(K,1); end
if numel(nk) == 1; nk = nk*ones(K,1); end
n = sum(nk);

[r, randomflag, orthflag, distribution] = getopts(options);

% generate labels
labels = zeros(1, n); cumsum_nk = 0;
for k = 1 : K
    labels(cumsum_nk+1 : cumsum_nk+nk(k)) = k;
    cumsum_nk = cumsum_nk+nk(k);
end


% generate data
X = zeros(m, n);

    
if r > sum(dk) || r <= max(dk); error('inappropriate r \n'); end
U = 1 - 2*rand(m); U = U(:, 1:r);
if orthflag; U = orth(U); end
P = zeros(r, sum(dk)); cumsum_dk = 0;
yk = cell(K, 1); uk = yk;
while 1
    for k = 1 : K
        Pk = 1-2*rand(r,dk(k));
        Uk = orth(U*Pk);
        Ik = labels == k;
        switch distribution
            case 'uniform'
                Yk = 1 - 2*rand(dk(k), nk(k));
            case 'normal'
                Yk = randn(dk(k), nk(k));
        end
        uk{k}=Uk; yk{k} = Yk;
        X(:, Ik) = Uk*Yk;
        % to calculate actual r, considering the possibility there is
        % a column of U not be selected
        P(:,cumsum_dk+1:cumsum_dk+dk(k)) = Pk;
        cumsum_dk = cumsum_dk + dk(k);
    end
    r0 = rank(P);
    if r0==r; break; end
end
    
    
    
% permute data
randidx = 1:n;
if randomflag
    randidx = randperm(n);
    X = X(:, randidx);
    labels = labels(randidx);
end

% generate settings
settings.randidx = randidx;

settings.r = r;
settings.uk = uk;
settings.Yk = yk;
settings.dk = dk;
settings.U = U;

% save data
if nargin == 6
    save(fname, 'X', 'labels', 'settings');
end


    function [r, randomflag, orthflag, distribution] = getopts(options)
        r = 1; randomflag = 0; orthflag = 1; distribution = 'uniform';
        if isfield(options,'r'); r = options.r; end
        if isfield(options, 'randomflag'); randomflag = options.randomflag; end
        if isfield(options, 'orthflag'); orthflag = options.orthflag; end
        if isfield(options, 'distribution'); distribution = options.distribution; end
    end
end
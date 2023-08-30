function prepare_create_data(foldname, REP)
% script to generate synthetic data in Section 8.2
% foldname: where we save these data
% REP: number of repeat of the experiment

% 10/04/2019 by Xia Yuqing


if ~exist('REP', 'var'); REP = 10; end
if ~exist(foldname, 'dir'); mkdir(foldname); addpath(foldname); end

% parameters
distribution = 'uniform';
Ranks = [10 14 20];
dcs = cell(max(Ranks),1); 
dcs([10 14 20]) = {4:8, 7:11, 11:15}; 
num_paras = 3*5;
paras = cell(num_paras,1);


i = 1;
for r = Ranks
    for dk = dcs{r}
        paras{i}.dk = dk;
        paras{i}.others.r = r;
        paras{i}.others.distribution = distribution;
        paras{i}.K = 5;
        paras{i}.nk = 50;
        paras{i}.m = 100;
        i = i+1;
    end
end

seedNo = 1;
fprintf('create data set with %s random seed = %d \n', distribution, seedNo);
rng(seedNo); 

for rep = 1 : REP
    for i = 1 : num_paras
        fname = ['paras', num2str(i), 'rep', num2str(rep),'.mat'];
        p = paras{i};
        [X, labels, settings] = gen_data(p.dk, p.nk, p.m, p.K, p.others);
        [idealC, idealV, idealW, V, V_bot] = getideal(X, settings.r, labels, settings.dk);
        save(fullfile(foldname, fname), 'X', 'labels', 'settings', 'idealC', 'idealV', 'idealW','V', 'V_bot');
    end
end

%fname = 'parasirepj.mat';
save(fullfile(foldname, 'syndata_log.mat'), 'REP', 'num_paras', 'fname', 'paras');

end % function

function [idealC, idealV, idealW, V, V_bot] = getideal(X, r, labels, dk)
[~, ~, V] = svd(X);
V_bot = V(:, r+1:end); V = V(:, 1:r);

K = numel(dk);
idealC = zeros(numel(labels));
idealV = zeros(numel(labels), sum(dk));
cumsum_dk = 0;
Vk = cell(K,1);
for k = 1 : K
    Ik = labels == k;
    Jk = cumsum_dk+1 : cumsum_dk+dk(k);
    cumsum_dk = cumsum_dk + dk(k);
    [~, ~, v] = svd(X(:, Ik), 'econ');
    Vk{k} = v(:, 1:dk(k));
    idealC(Ik, Ik) = Vk{k}*Vk{k}';
    idealV(Ik, Jk) = Vk{k};
end
% Calculate ideal W
F = idealV'*V; F_bot = null(F');
idealW = V_bot'*(idealV*F_bot);
end

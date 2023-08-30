% function gini_relv_gap
% The results will be managed by table3.m to get table3 in Section 8.2

% 10/04/2019 by Xia Yuqing

function gini_relv_gap(datafold, outfold)
% fold names
%datafold = 'syndata/gdata';
%outfold = 'out/gdata/';
if ~exist(outfold, 'dir'); mkdir(outfold); addpath(outfold); end

% algorithm names
num_algo = 8; 
algo = cell(num_algo,1);
algo{1} = 'MSS_MCG';
algo{2} = 'MSS_AO';
algo{3} = 'MSS_HO';
algo{4} = 'LRR';
algo{5} = 'SSC';
algo{6} = 'HardS3C';
algo{7} = 'SoftS3C';
algo{8} = 'LRSSC';

            
% initialize variants
load(fullfile(datafold, 'syndata_log.mat'));
conn = -ones(num_algo, num_paras, REP);
devi = -ones(num_algo, num_paras, REP);
rel_gap = -ones(num_algo, num_paras, REP);

for rep = 1 : REP
    fprintf('---------------REP = %d-------------------\n', rep);
    for i = 1 : num_paras
        fname = ['paras', num2str(i), 'rep', num2str(rep),'.mat'];
        load(fullfile(datafold, fname))
        K = numel(settings.dk);
        
        for j = 1 : num_algo  
            current_algo = algo{j};
            gname = fullfile(current_algo, fname);
            load(fullfile(outfold,gname));
            
            [~,~,out] = SpectralClustering(abs(C)+abs(C'),K);
            conn(j, i, rep) = IntraBConn(C, labels);
            devi(j, i, rep) = BdiagDevi(C, labels);
            rel_gap(j, i, rep) = out.rel_gap;
           
        end
    end
end

save(fullfile(outfold, 'gini_relv_gap.mat'), 'conn', 'devi', 'rel_gap');
            
            
            
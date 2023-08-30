% function test_gdata
% calculate C and missrate for data in the datafold via different
% algorithms, outputs are saved in outfold.

% The results will be managed by table3.m to get table3 in Section 8.2

% 10/04/2019 by Xia Yuqing

clear; %clc;

% fold names
datafold = 'syndata/gdata';
outfold = 'out/gdata/';
if ~exist(datafold, 'dir'); mkdir(datafold); addpath(datafold); end
if ~exist(outfold, 'dir'); mkdir(outfold); addpath(outfold); end

prepare_create_data(datafold)

saveflag = 1;

% algorithm names
num_algo = 6; 
algo = cell(num_algo,1);
algo{1} = 'MSS_MCG';
algo{2} = 'MSS_AO';
algo{3} = 'MSS_HO';
algo{4} = 'LRR';
algo{5} = 'SSC';
algo{6} = 'LRSSC';

s = ['mkdir ', outfold];
system(s);
for i = 1 : num_algo
    fdir = fullfile(outfold, algo{i});
    s = ['mkdir ', fdir];
    system(s);
end

% paras for MCG
MCGopts = [];
% paras for omega_mcg
AOopts.opts = MCGopts;
AOopts.lambda = 50;
% paras for LRSSC
LRSCopts.lambda = 50;
LRSCopts.ADAPTIVE_ON = 1; 
% paras for SSC
SSCopts.affine = 0;
SSCopts.alpha = 400;
% paras for strSSC
S3Copts.affine = 0;
S3Copts.outliers = 0;
S3Copts.T = 1;
S3Copts.iter_max = 10;
S3Copts.nu = 1;
S3Copts.r = 0;
S3Copts.lambda = 800;
S3Copts.mu_max = 1e8;
S3Copts.tol = 1e-3; % 1e-5
S3Copts.maxIter = 150;
S3Copts.rho = 1.1;
S3Copts.gamma0 = 0.1;
S3Copts.SSCrho = 1;
S3Copts.DEBUG = 0;

% S3Copts.gamma0 = 0.1;
% S3Copts.epsilon = 1e-3; % 1e-3

HardS3Copts = S3Copts;
HardS3Copts.sc_method = 'StrSSC-Fix2+';
HardS3Copts.gamma0 = 0.1;

SoftS3Copts = S3Copts;
SoftS3Copts.sc_method = 'softStrSSC-Fix2+';
SoftS3Copts.gamma0 = 0.1;
% paras for combine
HOopts.flag_display = 0;
            
% initialize variants
load(fullfile(datafold, 'syndata_log.mat'));
all_missrate = -ones(num_algo, num_paras, REP);
ts = -ones(num_algo, num_paras, REP);

numIter = zeros(3, num_paras, REP);
       
%%
for rep = 1 : REP
    fprintf('---------------REP = %d-------------------\n', rep);
    for i = 1 : num_paras
        fname = ['paras', num2str(i), 'rep', num2str(rep),'.mat'];
        load(fullfile(datafold, fname))
        lambda = 50;
        d = sum(settings.dk);  K = numel(settings.dk);
        r = settings.r; n = size(X, 2);
        [~,~,V] = svd(X);
        V_bot = V(:, r+1:end);
        V = V(:,1:r);
        
        Cv = V*V';
        W0 = [eye(d-r); zeros(n-d, d-r)];
        Omega = ones(n)-eye(n);
        
           
        for j = 1 : num_algo
            tic
            current_algo = algo{j};
            switch current_algo
                case 'MSS_MCG'
                    W = MSS_MCG(Omega,V,V_bot,lambda, W0, MCGopts);
                    VW = V_bot*W;
                    C = Cv + VW*VW';
                    [grps,~,~] = SpectralClustering(abs(C)+abs(C)', K);
                 case 'MSS_AO'
                      [grps, C, W , out] = MSS_AO(V,V_bot, d, K, AOopts);
                 case 'HardS3C'
                      [acc,  ~, C, ~] = StrSSCplus(X, labels, HardS3Copts);
                      all_missrate(j,i,rep) = 1-acc;
                 case 'SoftS3C'
                      [acc,  ~, C, ~] = StrSSCplus(X, labels, SoftS3Copts);
                      all_missrate(j,i,rep) = 1-acc;
                 case 'strSSC'
                      [acc,  ~, C, ~] = StrSSC(X, labels, S3Copts);
                      all_missrate(j,i,rep) = 1-acc;
                 case 'SSC'
                    C = admmLasso_mat_func(X, SSCopts.affine, SSCopts.alpha);
                    [grps,~,~] = SpectralClustering(abs(C)+abs(C)', K);
                case 'LRR'
                   C = Cv;
                   [grps,~,~] = SpectralClustering(abs(C)+abs(C)', K);
                case 'LRSSC'
                   C = ALM_LRSSC(X, LRSCopts.lambda, LRSCopts.ADAPTIVE_ON);
                   [grps,~,~] = SpectralClustering(abs(C)+abs(C)', K);
                case 'MSS_HO'
                   [grps, C, out] = MSS_HO(X, V, V_bot, d, K, lambda, HOopts,labels);
            end
            
            switch current_algo
                case {'MSS_AO', 'MSS_HO'}
                    numIter(j,i,rep) = out.iter;
                    all_missrate(j,i,rep) = Misclassification(grps, labels);
                case {'HardS3C','SoftS3C','strSSC'}
                    fprintf('');
                otherwise
                    all_missrate(j,i,rep) = Misclassification(grps, labels);
            end
                    
           
            ts(j ,i, rep) = toc;
            
             fprintf(' done for algo: %s params %d \n', current_algo, i );
             gname = fullfile(current_algo, fname);
             if saveflag; save(fullfile(outfold,gname), 'C'); end
            
        end
    end
end

if saveflag
    save(fullfile(outfold, 'missrate.mat'), 'all_missrate', 'ts',...
        'MCGopts', 'AOopts','HOopts', 'LRSCopts', 'HardS3Copts','SoftS3Copts','numIter');
end


            
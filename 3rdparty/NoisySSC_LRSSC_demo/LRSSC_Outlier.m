%--------------------------------------------------------------------------
% Extended version of LRSSC, i.e.
%   min_C 1/(1+lambda)||C||_1 + lambda/(1+lambda) ||C||_* + gamma||Y-YC||_1
%   s.t.  diag(C) = 0
% This function is based on admmOutlier_mat_func.m
%
% By ADMM, we rewrite the objective function to be
%   min lambda1 ||C1||_1 + lambda2 ||C2||_* + gamma ||E||_1
%   s.t. Y-YZ-E=0(mu1, Lambda1), Z-C1=0(mu2, Lambda2), Z-C2=0(mu3, Lambda3)
%        diag(Z) = diag(C1) = diag(C2) = 0
%
% Xiayq @ 20190614
%--------------------------------------------------------------------------

function [C, E] = LRSSC_Outlier(Y,lambda,alpha,thr,maxIter)

if (nargin < 2)
    lambda = 0;
end
if (nargin < 3)
    % default regularizarion parameters
    alpha = 20;
end
if (nargin < 4)
    % default coefficient error threshold to stop ALM
    % default linear system error threshold to stop ALM
    thr = 2*10^-4; 
end
if (nargin < 5)
    % default maximum number of iterations of ALM
    maxIter = 150; 
end

if (length(alpha) == 1)
    alpha1 = alpha(1);
    alpha2 = alpha(1);
    alpha3 = alpha(1);
elseif (length(alpha) == 2)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    alpha3 = alpha(2);
elseif (length(alpha) == 3)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    alpha3 = alpha(3);
end

if (length(thr) == 1)
    thr1 = thr(1);
    thr2 = thr(1);
    thr3 = thr(1);
elseif (length(thr) == 2)
    thr1 = thr(1);
    thr2 = thr(2);
    thr3 = thr(2);
end

lambda1 = 1/(1+lambda);
lambda2 = lambda/(1+lambda);

[D,N] = size(Y);

gamma = alpha3 / norm(Y,1);
P = [Y/lambda1 eye(D)/gamma];

% setting penalty parameters for the ADMM
mu1 = alpha1 * 1/computeLambda_mat(Y,P);
mu2 = alpha2 * 1;
mu3 = mu2;


% initialization
A = inv(mu1*(P'*P)+(mu2+mu3)*eye(N+D));
C1 = zeros(N+D,N);
C2 = zeros(N+D,N);

Lambda1 = zeros(D,N);
Lambda2 = zeros(N+D,N);
Lambda3 = zeros(N+D,N);
err1 = 10*thr1; err2 = 10*thr2; err3 = 10*thr3;
i = 1;
% ADMM iterations
while ( (err1(i) > thr1 || err2(i) > thr2 || err3(i) > thr3) && i < maxIter )
    % updating Z
    Z = A * (mu1*P'*(Y+Lambda1/mu1)+mu2*(C1-Lambda2/mu2)+mu3*(C2-Lambda3/mu3));
    Z(1:N,:) = Z(1:N,:) - diag(diag(Z(1:N,:)));
    % updating C1
    C1 = soft_thresh(Z+Lambda2/mu2, lambda1/mu2);
    C1(1:N,:) = C1(1:N,:) - diag(diag(C1(1:N,:)));
    % updating C2
    C2(1:N,:) = sigma_soft_thresh(Z(1:N,:)+Lambda3(1:N,:)/mu3, lambda2/(mu3*lambda1));
    C2(N+1:end,:) = Z(N+1:end,:);
    % updating Lagrange multipliers
    Lambda1 = Lambda1 + mu1 * (Y - P * Z);
    Lambda2 = Lambda2 + mu2 * (Z - C1);
    Lambda3 = Lambda3 + mu3 * (Z - C2);
    % computing errors
    err1(i+1) = errorCoef(Z,C1);
    err2(i+1) = errorCoef(Z,C2);
    err3(i+1) = errorLinSys(P,Z);
    %
    i = i + 1;
end
%fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),i);
C = C1(1:N,:)/lambda1;
E = C1(N+1:end,:)/gamma;

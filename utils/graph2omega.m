function [Omega, grps] = graph2omega(S, K, tau)
% the input S is a graph, which is after post-processing

    n = size(S,1);
    sd = diag(1./sqrt(sum(S)+eps));
    Laplacian = eye(n) - sd*S*sd;
    [~,~,Y] = svd(Laplacian);
    Y = Y(:,n-K+1:n); 
 
    % kmeans clustering
    MAXiter = 1000; % Maximum number of iterations for KMeans
    REPlic = 20; % Number of replications for KMeans
    norm_Y = sqrt(sum(Y.^2,2));
    Y = bsxfun(@rdivide, Y, norm_Y+eps);
    [grps,~,~,D] = kmeans(Y,K,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    D = sqrt(D'); % D: squraed distences of points to every centriod
            
    % Update Omega
    Dmin = min(D, [], 1);
    Dmax = max(D, [], 1);
    q = bsxfun(@minus, D, Dmin);
    q = bsxfun(@rdivide, q, Dmax-Dmin);
            
           
    P = real(q<tau); % tau = 0.3+0.2; %tau = 0.3+0.2/rep;
    P = bsxfun(@rdivide, P, sum(P,1));
    Omega = P'*P<1; %Omega = 1-P'*P>0.2;
end

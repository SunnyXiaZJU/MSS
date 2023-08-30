function [J, Dist_S, out]=SubSpaceEsti2(X,d,K,J,maxiter,weight)
% J is a colomn vector 
% modified from zhenyue zhang's first version--SubSpaceEsti.m

n = size(X,2);
SVDs = cell(K,1); Dist_S = zeros(maxiter,n);
for iter = 1:maxiter
    
    % SVDs of X_k
    s_ind = zeros(1,n); ei = 0; 
    for k=1:K
        Jk = J==k;
        Xk = X(:,Jk);
        [Gk,Sk,Qk] = svd(Xk,0); sk = diag(Sk);
        SVDs{k} = {Gk,sk,Qk};
        si = ei+1; ei = ei+length(sk);
        s_ind(si:ei) = k;
    end
    s_ind(ei+1:end) = [];
    
    % Weights and the updatings
    omega = zeros(K,1);
    switch weight
        case 'mean error'
            for k=1:K
                omega(k) = 1/sum(J==k);
            end
        case 'relative F-error'
            for k=1:K
                sk = SVDs{k}{2};
                omega(k) = 1/sum(sk.^2);
            end                    
        case 'relative 2-error'
            for k=1:K
                sk = SVDs{k}{2};
                omega(k) = 1/sk(1)^2;
            end  
        otherwise
            omega = ones(K,1);
    end
    
    % Determine implicitly
    ei = 0; ws = zeros(n,1);
    for k=1:K
        sk = SVDs{k}{2};
        si = ei+1; ei = ei+length(sk);
        ws(si:ei) = omega(k).*sk;
    end
    ws(ei+1:end) = [];
    [~,ind] = sort(ws,'descend');
    s_ind = s_ind(ind(1:d));
    
%     if numel(unique(s_ind)) < K; 
%         break; 
%     end
%     
    % Update J
    J_old = J;
    err_s = zeros(K,1); dist_S = zeros(K,n); dk = -ones(K,1);
    for k=1:K
        dk(k) = sum(s_ind == k);
        Gk = SVDs{k}{1}(:,1:dk(k));
        err_s(k) = norm(SVDs{k}{2}(dk(k)+1:end),'fro');
        Ek = X-Gk*(Gk'*X);
        dist_S(k,:) = omega(k)*sum(Ek.^2,1);
    end
    [dist_s,J] = min(dist_S,[],1);
    J = J';
    Dist_S(iter,:) = dist_s;
    
    % Check convergence
    if sum(abs(J_old-J))==0
        break
    end
    
end
 
Dist_S(iter+1:end,:) = [];
out.dk = dk;
out.err_s = err_s;
out.dist_S = dist_S;


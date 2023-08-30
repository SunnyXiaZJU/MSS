function [C, out] = PRRSS(X, Omega, d, lambda, alpha, beta, type, opts)
% min ||Omega.*C||_1 + 0.5lambda||c||_2^2 + alpha ||X-XC||_ell +
% 0.5beta||C-QQ'||_F^2 s.t. Q \in S
% calculate Q,C alternatively with objective function decreasing
% for the subproblem of calculating C, we use ADMM, namely
% min ||Omega.*Z||_1 + 0.5lambda||z||_2^2 + alpha ||E||_ell +
% 0.5beta||C-QQ'||_F^2 s.t. C = Z, X = XC+E
% final version

[m, n] = size(X);

[outiter, initer, rho1, rho2,  gamma, epsilon1, epsilon2] = getopts(opts,X, beta);

% initialization
if ~isfield(opts, 'C'); C = randn(n); else C = opts.C ; end
Y1 = zeros(n); Y2 = zeros(m, n);
Q = [eye(d); zeros(n-d,d)];
XTX = X'*X;


fvals = -ones(outiter, 5);
for iter = 1 : outiter
    C_old = C;
    
    % update C (Z, E, C)
    for i = 1 : initer
        in_C_old = C;
        
        % update Z
        R = (ones(n)-eye(n))/rho1 + eye(n)/(rho1+lambda);
        Z = R .* wthresh(rho1*C+Y1, 's', Omega);
        
        %update E
        switch type
            case 'fro'
                E = (rho2*(X-X*C) - Y2)/(2*alpha + rho2);
            case '1'
                E = wthresh(X-X*C-Y2/rho2, 's', alpha/rho2);
            case '21'
                A = X-X*C-Y2/rho2;
                normA = sqrt(sum(A.*A));
                E = max(normA-alpha/rho2, 0) ./ normA;
                E = bsxfun(@times, A, E);
        end
        
        % update Y1, Y2, rho1, rho1
        Y1 = Y1 + rho1*(C-Z);
        Y2 = Y2 + rho2*(X*C+E-X);
        rho1 = rho1*gamma;
        rho2 = rho2*gamma;
        
        % update C
        C = (rho2*XTX+(rho1+beta)*eye(n))\(rho2*XTX - rho2*X'*E - X'*Y2 + rho1*Z - Y1 + Q*Q'*beta);
        if norm(in_C_old - C, 'fro') < epsilon2 * norm(in_C_old, 'fro')
            break;
        end
    end
    
     % update Q
    B = 0.5*(C+C');
    [evt, evl]  = eig(B);
    evl = diag(evl);
    dplus = min([sum(evl>0), d]);
    [evl, tid] = sort(evl, 'descend');
    evt = evt(:,tid(1:dplus));
    evl = evl(1:dplus);
    Q =  evt*diag(sqrt(evl));
   
    
    f = [norm(Omega(:).*C(:),1), norm(diag(C),2)^2, -1, norm(C-Q*Q','fro')^2];
    switch type
        case 'fro'
            f(3) = norm(X-X*C, 'fro')^2;
        case '1'
            f(3) = sum(sum(abs(X-X*C)));
        case '21'
            f(3) = sum(sqrt(sum(abs(X-X*C))));
    end
    fvals(iter, :) = [f(1)+0.5*lambda*f(2)+alpha*f(3)+0.5*beta*f(4), f];
    
    if iter > 3 && norm(C_old - C, 'fro') < epsilon1*norm(C_old,'fro')
        break;
    end
    
end
out.fvals = fvals(1:iter, :);

end


function [outiter, initer, rho1, rho2,  gamma, epsilon1, epsilon2] = getopts(opts, X, beta)
dd = svd(X, 'econ');
rho1 = beta / sqrt(dd(1)^2+1);
rho2 = rho1;
outiter = 100; 
gamma = 1.1; epsilon1 = 1e-5; initer = 5; epsilon2 = 1e-5;

fnames = fieldnames(opts);
for i = 1 : numel(fnames)
    if exist(fnames{i}, 'var')
        eval([fnames{i}, '=opts.', fnames{i}, ';']);
    end
end
end

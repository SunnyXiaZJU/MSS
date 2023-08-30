function [W,C, out] = MSS_MCG(Omega,V,V_bot,lambda,W,options)

% MSS_MCG algorithm for Section 8.2
% 10/04/2019 by Xia Yuqing

[delta, gamma, eps_d, maxiter_delta, maxiter_w, eps_w, alpha_init, rho, ...
    maxiter_alpha, eps_a, tau_armijo, criterion_W] = getopts(options);


C_v = V*V';
VW = V_bot*W;
C = C_v+VW*VW'; c = diag(C); 
fvals = -ones(maxiter_delta*maxiter_w,1);
count = 1;

%fprintf(' delta   ||P||_F   alpha  alpha*||H||_F iter_w mena_it_a  rho\n') 
for iter_delta=1:maxiter_delta
    alpha = alpha_init;
    C_old = C;
    
    % Initialization for minimize f_delta
    aOC = abs(Omega.*C); 
    fd = f_delta(aOC,c,delta,lambda);
    
    % Solve a minimizer of f_delta by MCG
    n_a = 0; 
    
    for iter_w=1:maxiter_w
        
        % Determe the manifold-restricted CG direction 
        sgn_C = sign(C);
        S = Omega.*sgn_C.*min(aOC/delta,1);
        %disp(norm(S,'fro'))
        gradf = 2*(V_bot'*(S+diag(lambda*c))*VW);
        P = proj_(W,gradf); 
        if iter_w == 1
            H = -P; 
        else
            Y = proj_(W,gradf_old); Y = gradf-Y;
            Z = proj_(W,H);
            b = sum(Y(:).*Z(:));
            if abs(b)>1.e-10
                Y = Y-(2*sum(Y(:).^2)/b)*Z;
                beta = sum(Y(:).*P(:))/b;
                H = -P+beta*Z;
            else
                disp('fast desending')
                H = -P;
            end
        end
        H = H/norm(H,'fro');
        
        % Update W via line search under Armijo condition
        p = tau_armijo*sum(P(:).*H(:));
        VH = V_bot*H;
                
        iter_alpha = 0;
        flag_1 = 0; flag_2 = 0; a0 = alpha;
        while iter_alpha<maxiter_alpha
            iter_alpha = iter_alpha+1;
            
            VW_new = VW+alpha*VH;
            C = C_v + VW_new*VW_new'; c = diag(C);
            aOC = abs(Omega.*C);
            fd_new = f_delta(aOC,c,delta,lambda);
            
            if fd_new <= fd+alpha*p
                flag_1 = 1;
                a_W = alpha;
                fd_W = fd_new;
                V_W = VW_new;
                C_W = C; c_W = c;
                aOC_W = aOC;
            else
                flag_2 = 1; %a_plus = alpha;
            end
            if flag_1*flag_2==1
                W = W+a_W*H; 
                fd_old = fd;
                fd = fd_W; VW = V_W; C = C_W; c = c_W; aOC = aOC_W;
%                f = sum(aOC(:))+sum(c.^2)*lambda/2;
%                Ra = [delta iter_w iter_alpha rho a0 a_W a_plus fd (1-f/fd)/delta];
%                str = ' delta= %7.7f it_w= %3d it_a= %2d rho= %4.1e a0= %4.1e a_W= %4.1e a_W+= %4.1e fd= %18.16f ref= %4.2f\n';
%                fprintf(str,Ra)
                alpha = a_W;          
                break
            elseif flag_2==0
                alpha = alpha/rho;
            else
                a_W = alpha;
                alpha = alpha*rho;
            end
            
        end
        n_a = n_a+iter_alpha;
        gradf_old = gradf;
        
        fvals(count) = fd;
        count = count + 1;
%        save(['zy/out/delta' num2str(delta) 'iterw' num2str(iter_w) '.mat'], 'C')
        % Check convergence of the inner iteration
        switch criterion_W
            case 'a_W'
                if a_W < eps_a, break, end
                
            case 'gap of f_delta'
                n_delta = sum(aOC(:)<=delta);
                if fd_old-fd <= 0.001*n_delta*delta, break, end
        end
    end
    
%    R = [delta iter_w n_a/iter_w];
%    fprintf(' delta= %7.7f it_w= %4d mean_it_a= %4.1f\n',R);
    
    % Check convergence of the outer iteration
    if max(abs(Omega(:).*(C(:)- C_old(:))))<eps_w && delta < eps_d
        out.fvals = fvals(1:count-1);
        return
    else
        delta = gamma*delta;
    end
    
end
out.fvals = fvals(1:count-1);
end

function f = f_delta(aOC,c,delta,lambda)
    
    a = aOC(:); idx = a>delta; 
    a1 = a(idx); a2 = a(~idx);
    f = sum(a1)+sum(a2.^2 + delta^2)/(2*delta) + sum(c.^2)*lambda/2;
    
end

function H = proj_(W,G)

    [Q,S,~] = svd(W'*W, 'econ'); s = diag(S);
    F = W'*G; F = Q'*(F-F')*Q;
    N = F./(bsxfun(@plus, s, s'));
    N = Q*N*Q';
    H = G-W*N; 
    
end

function [delta, gamma, eps_d, maxiter_d, maxiter_w, eps_w, alpha, ...
    rho, maxiter_a, eps_a, tau, crit] = getopts(opts)
   
    delta = 1; gamma = 0.1; eps_d = 1e-5; maxiter_d = 10;
    maxiter_w = 800; eps_w = 1e-5; alpha = 1;
    rho = 0.5; maxiter_a = 10; eps_a = 1e-3;
    tau = 0.5; crit = 'a_W';

%     delta = 1; gamma = 0.1; eps_d = 1e-5; maxiter_d = 10;
%     maxiter_w = 500; eps_w = 1e-4; alpha = 1;
%     rho = 0.5; maxiter_a = 10; eps_a = 1e-4;
%     tau = 0.5; crit = 'a_W';

    if ~isempty(opts)
        fnames = fieldnames(opts);
        for i = 1 : numel(fnames)
            if exist(fnames{i}, 'var')
                eval([fnames{i}, '=opts.', fnames{i}, ';']);
            end
        end
    end
            
       
end
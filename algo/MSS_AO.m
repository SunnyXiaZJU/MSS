function [grps, C, W , out] = MSS_AO(V,V_bot, d, K, options, labels)

% MSS_AO algorithm for Section 8.2
% 10/04/2019 by Xia Yuqing

out = [];
[n, r] = size(V);
[lambda0, Omega, outN, tau, flag_display, opts] = getopts(options, n);


DEBUG = 0;
if exist('labels', 'var'); DEBUG = 1; end

W0 =  [eye(d-r); zeros(n-d, d-r)];


 for iter = 1 : outN
            if iter == 1
                lambda = lambda0;
            else
                lambda = min(lambda0, 2*sum(abs( Omega(:).*C(:)))/norm(diag(C),2)^2);
            end
            Omega_old = Omega;
            
            % MCG solution
             
            [W,C] = MSS_MCG(Omega,V,V_bot,lambda, W0, opts);  
           
            % Spectral projection
            S = abs(C)+abs(C');
            sd = diag(1./sqrt(sum(S)+eps));
            Laplacian = sd*S*sd;
            [~,~,Y] = svd(Laplacian); 
            Y = Y(:,1:K); 
           
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
          
            P = real(q<tau);
            P = bsxfun(@rdivide, P, sum(P,1));
            Omega = sign(1-P'*P);
     
            % show Omega
            if flag_display
            figure(1),imshow(Omega);
            title(['prime: iter = ' num2str(iter)], 'FontSize', 20,'FontWeight','bold');
            pause(.1)
            end
            
            % classification error
            if DEBUG; mr = Misclassification(grps, labels); end
       
           
            if norm(Omega_old - Omega, 'fro') == 0; break; end
 end
 if DEBUG; out.mr = mr; end
 out.iter = iter;
 
end

function [lambda, Omega,outN, tau, flag_display, opts] = getopts(options, n)
    opts = []; lambda = 50; outN = 10; flag_display = 0; tau = 0.5;Omega = ones(n)-eye(n);
    if ~isempty(options)
         fnames = fieldnames(options);
         for i = 1 : numel(fnames)
            if exist(fnames{i}, 'var')
                eval([fnames{i}, '=options.', fnames{i}, ';']);
            end
        end
    end
end
        
 
 
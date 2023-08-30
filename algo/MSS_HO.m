function [grps1, C, out] = MSS_HO(X, V, V_bot, d, K, lambda0, options, labels)

% MSS_HO algorithm for Section 8.2
% 10/04/2019 by Xia Yuqing
    
    out = []; DEBUG = 0;
    if exist('labels', 'var'); DEBUG = 1; end
    [outN, maxiter_mcg, maxiter_ks, beta, tau, flag_display] = getopts(options);

    n = size(V,1);
  
    opts.outN = maxiter_mcg;
    opts.tau = tau;
    opts.flag_display = flag_display;
    opts.lambda = lambda0;
    
    Omega = ones(n)-eye(n);
 
    for outiter = 1 : outN
        Omega_old = Omega;
       
        opts.Omega = Omega;
        if outiter > 1
            opts.lambda = min(lambda0, 2*sum(abs( Omega(:).*C(:)))/norm(diag(C),2)^2);
        end
       
        if DEBUG
            [grps, C, ~, temp_out] = MSS_AO(V,V_bot, d, K, opts, labels);
            mr = temp_out.mr;
        else
            [grps, C, ~ , ~] = MSS_AO(V,V_bot, d, K, opts);
        end
     
        %%
        
        [grps1, ~, ~]=SubSpaceEsti2(X,d,K,grps,maxiter_ks, 'ones');
        if DEBUG
            mr1 = Misclassification(grps1, labels);
        end
        Omega = ones(n);
        for i = 1 : max(grps1)
            Ik = grps1 == i;
            Omega(Ik, Ik) = 0;
        end
        if  numel(unique(grps1)) < K 
            Omega = real(Omega);
            Omega(Omega == 0) =  beta;%0.8;
            Omega = Omega - diag(diag(Omega));
        end
         if flag_display
         figure(1), imshow(Omega, 'Border','loose'); 
         title(['dual' num2str(outiter)], 'FontSize', 20,'FontWeight','bold');pause(.1);
        end
        if norm(Omega_old-Omega,'fro') == 0; out.grps = grps; break; end
       
    end
    if DEBUG
        out.mr = mr;
        out.mr1 = mr1;
    end
    out.grps = grps;
    out.iter = outiter;
end



function [outN, maxiter_mcg, maxiter_ks, beta, tau, flag_display] = getopts(options)
    outN = 10; maxiter_mcg = 10; maxiter_ks = 10;
    beta = 0.8; tau = 0.5; flag_display = 1;
    
    if ~isempty(options)
        fnames = fieldnames(options);
        for i = 1 : numel(fnames)
            if exist(fnames{i}, 'var')
                eval([fnames{i}, '=options.', fnames{i}, ';']);
            end
        end
    end

end


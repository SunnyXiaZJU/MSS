function [Omega, gs, D] = getomega(C, K, delta, mode)
% delta gap
% mode   -- 'overlap'  free the edge that the two nodes has overlap
%        -- 'disconnect' 
[gs, D] = SpectralClustering(abs(C)+abs(C'), K);
delta = (max(D,[],2)-min(D,[],2))*delta + min(D,[],2);
grps = D < delta*ones(1,K);
grps = diag(1./sum(grps,2))*grps;
n = size(C,1);

% R = min(D,[],2)./max(D,[],2);
% disp([min(R),max(R)]);
% disp(size(R));
% return;

switch mode,
    case {'overlap', 0}
        Omega = grps*grps' == 0;
    case {'disconnect', 1}
        Omega = grps*grps' <1;
    case 'probability'
        Omega = zeros(n);
        for i = 1 : n
            for j = 1 : n
                %if gs(i) ~= gs(j) 
                    Omega(i,j) = exp(-D(i, gs(i))/D(i,gs(j)));
                %end
            end
        end
        %Omega = 0.5*( Omega + Omega');
        Omega = Omega.*Omega';
    case 'disover'
        Omega1 = grps*grps' == 0;
        Omega2 = grps*grps' < 1;
        Omega = 0.5*Omega1+0.5*Omega2;
        
    otherwise
        Omega = ones(n);
        for i = 1 : K
            Ik = gs == i;
            Omega(Ik,Ik) = 0;
        end
end

end


function [groups,D] = SpectralClustering(CKSym,n)
% warning off;
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans
REPlic = 20; % Number of replications for KMeans

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}

DN = diag( 1./sqrt(sum(CKSym)+eps) );
LapN = speye(N) - DN * CKSym * DN;
[~, ~, vN] = svd(LapN);
kerN = vN(:,N-n+1:N);
for i = 1:N
    kerN(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
end
[groups, ~,~,D] = kmeans(kerN,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
end



%This code implements face clustering experiment in IEEE signal processing
%letter paper 'Robust Subspace Clustering via Smoothed Rank Approximation'
% zhao kang, 08/2015, zhao.kang@siu.edu


load YaleBCrop025.mat;


rho=.0003; % two tuning parameters
betazero=.4;

alpha=2;  % W construct
type=1;% different modeling of E
iter=100;

nSet = [2 3 5 8 10];
for i = 1 : length(nSet)
    n = nSet(i);
    idx = Ind{n};
    for j = 1 : size(idx,1)
        X = [];
        for p = 1 : n
            X = [X Y(:,:,idx(j,p))];
        end
        X = mat2gray(X);
        [m,nn]=size(X);
        
        %Initializations
        W=eye(nn);
        Y1=zeros(m,nn);
        Y2=zeros(nn);
        E=zeros(m,nn);
        
        beta=betazero;
        
        for ii=1:iter
            if ii>=2
                Zold=Z;
            end
            
            gamma=W-Y2/beta;
            [J]=LogSquare(gamma,beta/2);
            [E]=errormin(Y1,X,W,rho,beta,type);
            W=inv(eye(size(J))+X'*X)*(X'*(X-E)+J+(X'*Y1+Y2)/beta);
            Y1=Y1+beta*(X-X*W-E);
            Y2=Y2+beta*(J-W);
            beta=1.1*beta;
            Z=W;
            
            if type==1
                enorm=sum(sum(abs(X-X*Z)));
            elseif type==2
                enorm=sum(sum((X-X*Z).^2));
            else
                enorm=sum(sqrt(sum((X-X*Z).^2,1)));
            end
            
            if ii>3 && norm(Z-Zold,'fro')/norm(Zold,'fro')<1e-5
                break
            end
            
            func(ii)=log(det(eye(nn)+Z'*Z))+rho*enorm;
        end
        
        
        [U ss V] = svd(Z);
        ss = diag(ss);
        r = sum(ss>1e-6);
        
        U = U(:, 1 : r);
        ss = diag(ss(1 : r));
        V = V(:, 1 : r);
        
        M = U * ss.^(1/2);
        mm = normr(M);
        rs = mm * mm';
        L = rs.^(2 * alpha);
        
        actual_ids = spectral_clustering(L, n);
        err = calculate_err(s{n}, actual_ids);
        disp(err);
        
        
        missrateTot{n}(j) = err;
        
        %     dlmwrite('faceclustering.txt', [n rho beta alpha err] , '-append', 'delimiter', '\t', 'newline', 'pc');
        
    end
    avgmissrate(n) = mean(missrateTot{n});
    medmissrate(n) = median(missrateTot{n});
    %dlmwrite('face.txt', [n rho betazero alpha avgmissrate(n) medmissrate(n)] , '-append', 'delimiter', '\t', 'newline', 'pc');
    disp([n rho betazero alpha avgmissrate(n) medmissrate(n)]);
    save('out/faces/logdet.mat','missrateTot', 'avgmissrate', 'medmissrate');
end
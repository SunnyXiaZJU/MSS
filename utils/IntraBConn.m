function val = IntraBconn(C, labels)
% IntarBConn defined in (60) as 
% \sum_{\ell =1}^M \frac{|\tilde c_\ell}{\|C_{J^*}\|_1}(\frac{2(M-ell)+1}{M})

% It is between [0,1] equal to 1-GiniIndex. 
% It goes larger when the connection is better

% 11/04/2019


M = 0;
for k = unique(labels),
    M = M + sum(labels==k)^2;
end

i = 0;
c = ones(M,1);
for k = unique(labels),
    Ik = labels==k;
    j = sum(Ik)^2 ; 
    c (i +(1:j) ) = C(Ik,Ik);
    i = i+j;
end

c = sort(abs(c));
val = 2* sum(c/ sum(c) .* ((M -(1:M)'+ 0.5 )/M));
end


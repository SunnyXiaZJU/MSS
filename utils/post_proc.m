function G = post_proc(Z, beta, tau, alpha)
% construct graph by representatin matrix Z with parameter
% beta : sparsity-truncated parameter
% sigma: rank-truncated parameter
% alpha : sigmoid power


 Z = BuildAdjacency(thrC(Z,beta));
    
  [U, s, ~] = svd(Z);
  s = diag(s);
  r = sum(s>tau);
    
  U = U(:, 1 : r);
  s = diag(s(1 : r));
    
    
   M = U * s.^(1/2);
   mm = normr(M);
   rs = mm * mm';
   G = rs.^alpha;
    
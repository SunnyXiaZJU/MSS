function Y=sigma_soft_thresh(X,beta)
[U, S, V]=svd(X,'econ');
s = diag(S);
idx = find(s>beta);
U = U(:,idx); V = V(:, idx);
s = s(idx)-beta;
Y = U*diag(s)*V';
end
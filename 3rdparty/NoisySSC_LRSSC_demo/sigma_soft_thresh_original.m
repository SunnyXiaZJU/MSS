function Y=sigma_soft_thresh_original(X,beta)
[U S V]=svd(X,'econ');
v=soft_thresh_original(diag(S),beta);
Y=U*diag(v)*V';
end
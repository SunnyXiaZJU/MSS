function [U, iter] = svdr(Z, U, maxIter, tol)
if ~exist('maxIter', 'var'), maxIter = 50; end
if ~exist('tol', 'var'), tol = 1e-4; end
 for iter = 1: maxIter
    U_old = U;
    [U, ~] = qr(Z*U,0);
    ss = svd(U'*U_old);
    err_u = sqrt(1-min(ss)^2);
    if err_u<tol, break; end
end


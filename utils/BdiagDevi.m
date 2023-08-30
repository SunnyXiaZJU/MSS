function rel = BdiagDevi( C, labels )
% BdiagDevi defined in (59) as
% \|C_off(J)\|_1 / (\|C\|_1 - \|diag(C)\|_1)

% 11/04/2019

C = abs(C);
val = 0;

for k = unique(labels)
    Ik = labels==k;
    val= val + sum(sum(C(Ik, Ik)));
end

normC = sum(C(:));
den = normC - sum(diag(C));
if den == 0; den = den+eps; end
rel = (normC - val) /den ;

end


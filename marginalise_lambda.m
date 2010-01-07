function [rhos, rs, lambdas, pxd, p] = marginalise_lambda(file, olambda, op)

% file should contain lambdas, rhos, rs and pxd
load(file);

if nargin <= 1
    p = squeeze(trapz(rhos, trapz(rs, pxd)));
else
    p = interp1(olambda, op, lambdas, 'cubic', 0.0);
end

pxd = trapz(lambdas, pxd .* repmat(reshape(p, [1 1 numel(lambdas)]), numel(rs), numel(rhos)), 3);

end
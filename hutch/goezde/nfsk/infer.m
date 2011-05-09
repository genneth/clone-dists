function infer

mcm = [0 0 24 7 6 4 3 0 1 1 1 0 0 0 0 1 1 0 2 0 1 1 1 2 1 1 1 0 0 0 0 1 0 0 1 1 0 0 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
% ki67 = [0 0 26 7 8 4 3 1 3 2 2 3 0 1 1 2 3 2 0 1 1 1 0 0 3 0 1 0 0 0 1];
% edu = [0 0 50 14 14 8 6 1 4 3 3 3 0 1 1 3 4 2 2 1 2 2 1 2 4 1 2 0 0 0 1 1];

% use the mcm data to generate the parameters
rs = linspace(0.01, 0.49, 100);
lps = zeros(size(rs));
for i=1:numel(rs)
    theory = terminal_clone_dist_noloss(rs(i),numel(mcm)-1);
    lps(i) = dot(mcm(3:end), log(theory(3:end)));
end
ps = exp(lps - max(lps));
plot(rs, ps);
Z = sum(ps);
fprintf('r = %.2f +- %.2f\n', dot(rs, ps)/Z, dot(rs.*rs, ps)/Z);

end
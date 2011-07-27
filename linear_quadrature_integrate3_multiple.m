function results = linear_quadrature_integrate3_multiple(integrands, samples, pf)

% we assume that samples is a cell array, with structure:
% samples{i,:} = {r, gamma, lambda, [log p2], [log p3]}
% integrand :: (r, gamma, lambda) -> Double

coords = [[samples{:,1}]' [samples{:,2}]' [samples{:,3}]'];
warning off MATLAB:delaunay3:DeprecatedFunction;
tetras = delaunay3(coords(:,1), coords(:,2), coords(:,3));
lps = cellfun(pf, samples(:,4), samples(:,5)) - cellfun(@(l) log(l(1)), samples(:,3));
ps = exp(lps - max(lps)); % normalised to max(ps) == 1

v1 = coords(tetras(:,1),:);
v2 = coords(tetras(:,2),:); dv2 = v2 - v1;
v3 = coords(tetras(:,3),:); dv3 = v3 - v1;
v4 = coords(tetras(:,4),:); dv4 = v4 - v1;

volumes = 1/6 * abs(dot(dv2, cross(dv3,dv4, 2), 2));

vv = cellfun(@(i) feval(i, coords(:,1), coords(:,2), coords(:,3)), integrands, 'UniformOutput', false);

results = zeros(size(integrands));

for i = 1:numel(integrands)
    v = vv{i};
    results(i) = 1/24 * sum(sum(v(tetras) .* ps(tetras), 2) ./ volumes);
end

end

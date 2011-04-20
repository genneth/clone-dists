function result = linear_quadrature_integrate3(integrand, samples, pf)

% we assume that samples is a cell array, with structure:
% samples{i,:} = {r, gamma, lambda, [log p2], [log p3]}
% integrand :: (r, gamma, lambda) -> Double

coords = [[samples{:,1}]' [samples{:,2}]' [samples{:,3}]'];
warning off MATLAB:delaunay3:DeprecatedFunction;
tetras = delaunay3(coords(:,1), coords(:,2), coords(:,3));
lps = cellfun(pf, samples(:,4), samples(:,5)) - cellfun(@(l) log(l(1)), samples(:,3));
ps = exp(lps - max(lps)); % normalised to max(ps) == 1

[ntetras, ~] = size(tetras);
result = 0;

for i = 1:ntetras
    vertices = coords(tetras(i,:),:);
    y = zeros(4,1);
    for j = 1:4
        y(j) = ps(tetras(i,j)) * ...
            integrand(...
                coords(tetras(i,j),1), ...
                coords(tetras(i,j),2), ...
                coords(tetras(i,j),3));
    end
    volume = 1/6 * abs(det([
        vertices(2,:) - vertices(1,:);
        vertices(3,:) - vertices(1,:);
        vertices(4,:) - vertices(1,:)]));
    result = result + 1/24*sum(y)/volume;
end

end

function result = linear_quadrature_integrate3(integrand, samples)

% we assume that samples is a cell array, with structure:
% samples{i,:} = {prob, r, gamma, lambda, [log p2], [log p3]}
% integrand :: (r, gamma, lambda) -> Double

coords = [[samples{:,2}]' [samples{:,3}]' [samples{:,4}]'];
tetras = delaunay3(coords(:,1), coords(:,2), coords(:,3));
lps = cellfun(@(p2s, p3s)(sum(p2s) + sum(p3s)), {samples{:,5}}, {samples{:,6}});
ps = exp(lps - max(lps)); % normalised to max(ps) == 1

[ntetras, ~] = size(tetras);
result = 0;

for i = 1:ntetras
    vertices = zeros(4,3);
    y = zeros(4,1);
    for j = 1:4
        vertices(j,:) = coords(tetras(i,j),:);
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

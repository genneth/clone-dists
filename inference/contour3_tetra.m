function contour3_tetra(ah, X, Y, Z, W, c)

% we assume that samples is a cell array, with structure:
% samples{i,:} = {prob, r, gamma, lambda, [log p2], [log p3]}
% integrand :: (r, gamma, lambda) -> Double

tetras = delaunay3(X, Y, Z);

[ntetras, ~] = size(tetras);

triangles = zeros(0,3);
tverts = zeros(0,3);
n = 0;
for i = 1:ntetras
    vertices = [X(tetras(i,:)) Y(tetras(i,:)) Z(tetras(i,:))];
    y = W(tetras(i,:));
    if min(y) > c || max(y) < c
        continue; % contour does not go through this tetrahedron
    end
    [y, SI] = sort(y);
    vertices = vertices(SI,:);
    if y(1) > c || y(4) < c
        continue;
    elseif y(2) > c && y(3) > c
        tverts(end+1,:) = (vertices(2,:) * (c - y(1)) + vertices(1,:) * (y(2) - c)) / (y(2) - y(1));
        tverts(end+1,:) = (vertices(3,:) * (c - y(1)) + vertices(1,:) * (y(3) - c)) / (y(3) - y(1));
        tverts(end+1,:) = (vertices(4,:) * (c - y(1)) + vertices(1,:) * (y(4) - c)) / (y(4) - y(1));
        triangles(end+1,:) = [n+1 n+2 n+3];
        n = n+3;
    elseif y(2) < c && y(3) < c
        tverts(end+1,:) = (vertices(1,:) * (c - y(4)) + vertices(4,:) * (y(1) - c)) / (y(1) - y(4));
        tverts(end+1,:) = (vertices(2,:) * (c - y(4)) + vertices(4,:) * (y(2) - c)) / (y(2) - y(4));
        tverts(end+1,:) = (vertices(3,:) * (c - y(4)) + vertices(4,:) * (y(3) - c)) / (y(3) - y(4));
        triangles(end+1,:) = [n+1 n+2 n+3];
        n = n+3;
    else
        tverts(end+1,:) = (vertices(3,:) * (c - y(1)) + vertices(1,:) * (y(3) - c)) / (y(3) - y(1));
        tverts(end+1,:) = (vertices(4,:) * (c - y(1)) + vertices(1,:) * (y(4) - c)) / (y(4) - y(1));
        tverts(end+1,:) = (vertices(2,:) * (c - y(3)) + vertices(3,:) * (y(2) - c)) / (y(2) - y(3));
        tverts(end+1,:) = (vertices(2,:) * (c - y(4)) + vertices(4,:) * (y(2) - c)) / (y(2) - y(4));
        if dot(cross(tverts(n+2,:)-tverts(n+1,:),tverts(n+3,:)-tverts(n+1,:)), ...
                cross(tverts(n+3,:)-tverts(n+1,:),tverts(n+4,:)-tverts(n+1,:))) > 0
            triangles(end+1,:) = [n+1 n+2 n+3];
            triangles(end+1,:) = [n+1 n+3 n+4];
        elseif dot(cross(tverts(n+2,:)-tverts(n+1,:),tverts(n+4,:)-tverts(n+1,:)), ...
                cross(tverts(n+4,:)-tverts(n+1,:),tverts(n+3,:)-tverts(n+1,:))) > 0
            triangles(end+1,:) = [n+1 n+2 n+4];
            triangles(end+1,:) = [n+1 n+4 n+3];
        else
            triangles(end+1,:) = [n+1 n+2 n+3];
            triangles(end+1,:) = [n+1 n+4 n+2];
        end
        n = n+4;
    end
end

axes(ah);
trisurf(triangles, tverts(:,1), tverts(:,2), tverts(:,3), ...
    'FaceVertexAlphaData', 0.6 * ones(n,1), 'FaceAlpha', 'flat');

end

function Pb = Pb(P)

sz = size(P);
Pb = zeros(sum(sz)-1, 1);

for i = [1:sz(1)]
    for j = [1:sz(2)]
        Pb(i+j-1) = Pb(i+j-1) + P(i,j);
    end
end

end
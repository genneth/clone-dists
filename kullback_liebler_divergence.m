% compute D_KL(P || Q) = sum P_i log(P_i/Q_i)
function D = kullback_liebler_divergence(P, Q)

assert(all(size(P) == size(Q)));
% assert(P > 0 & Q > 0);

D = 0;

for i = 1:numel(P)
    if P(i) > 0
        D = D + P(i) * log(P(i) / Q(i));
    end
end

end
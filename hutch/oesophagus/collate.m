function c = collate(vs)
    m = max(cellfun(@numel, vs));
    c = zeros(1,m);
    for a=1:numel(vs)
        c(1:numel(vs{a})) = c(1:numel(vs{a})) + vs{a};
    end
end


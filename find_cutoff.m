function cutoff = find_cutoff(pxd, p)

ps = reshape(pxd, 1, numel(pxd));
cps = cumsum(sort(ps, 1, 'descend'));
n = find(cps > p*cps(end), 1);
% the right cutoff is now between n-1 and n'th elements

cutoff = mean([ps(n-1) ps(n)]);

end
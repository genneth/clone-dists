function cutoff = find_cutoff(pxd, p)

ps = sort(reshape(pxd, numel(pxd), 1),1,'descend');
cps = cumsum(ps);
n = find(cps > p*cps(end), 1);
% the right cutoff is now between n-1 and n'th elements

cutoff = mean([ps(n-1) ps(n)]);

end
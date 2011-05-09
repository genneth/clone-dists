function ps = terminal_clone_dist_noloss(r, N)

assert(N>2, 'must ask for more than the first 2 p(n)s');

ps = zeros(N+1,1);
ps(2+1) = r;
for i=3:N
    js = 2:(i-2);
    ps(i+1) = r*dot(ps(js+1), ps(fliplr(js)+1)) + (1-2*r)*ps(i);
end

end
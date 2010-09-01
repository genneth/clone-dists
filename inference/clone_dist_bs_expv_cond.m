function pbs = clone_dist_bs_expv_cond(Tl, Tr, Tg, P0, r, gamma, ts)

dist3_ = clone_dist_bs_expv(Tl, Tr, Tg, P0, r, gamma, ts);
basal_ = clone_dist_b(r, gamma, ts, 4);
[~,~,l_] = size(dist3_);
for j = 1:l_
    Z = 1 - basal_(0+1,j) - dist3_(1+1,0+1,j);
    dist3_(0+1,:,  j) = 0;
    dist3_(1+1,0+1,j) = 0;
    dist3_(:,:,j) = dist3_(:,:,j) ./ Z;        
end

pbs = dist3_;

end

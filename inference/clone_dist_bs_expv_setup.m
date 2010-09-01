function [Tlambda Tr Tgamma P0] = clone_dist_bs_expv_setup(k)

% The transition matrix only really depends on the max number of cells we
% will bother to track. This is good.

% k is the maximum number of cells we track

p = addpath(strcat(pwd, '/clone_dist_bs_expv_private'));

M = tetra(k+1);
Tlambda = sparse(M,M);
Tr      = sparse(M,M);
Tgamma  = sparse(M,M);

for i = 1:M
    t1 = i2t(i);
    
    % AA
    t2 = t1 - [1 0 0];
    if(valid_triple(t2))
        Tr(i, t2i(t2)) = t2(1);
    end
    
    % AB
    t2 = t1 - [0 1 0];
    if(valid_triple(t2))
        Tlambda(i, t2i(t2)) = t2(1);
        Tr(i, t2i(t2)) = -2 * t2(1);
    end
    
    % BB
    t2 = t1 - [-1 2 0];
    if(valid_triple(t2))
        Tr(i, t2i(t2)) = t2(1);
    end
    
    % C
    t2 = t1 - [0 -1 1];
    if(valid_triple(t2))
        Tgamma(i, t2i(t2)) = t2(2);
    end
    
    % exits
    Tlambda(i, i) = -t1(1);
    Tgamma(i, i)  = -t1(2);
    
end

P0 = zeros(M, 1);
P0(t2i([1 0 0])) = 1;

path(p);

end

function [Tlambda Tr Tgamma] = generate_transition_matrix(sz)

M = prod(sz);
Tlambda = sparse(M,M);
Tr      = sparse(M,M);
Tgamma  = sparse(M,M);

for i=[1:M]
    p1 = i2p(sz, i);
    
    % AA
    p2 = p1 - [1 0];
    if(valid_pair(sz, p2))
        Tr(i, p2i(sz, p2)) = p2(1);
    end
    
    % AB
    p2 = p1 - [0 1];
    if(valid_pair(sz, p2))
        Tlambda(i, p2i(sz, p2)) = p2(1);
        Tr(i, p2i(sz, p2)) = -2 * p2(1);
    end
    
    % BB
    p2 = p1 - [-1 2];
    if(valid_pair(sz, p2))
        Tr(i, p2i(sz, p2)) = p2(1);
    end
    
    % B -> nothing
    p2 = p1 - [0 -1];
    if(valid_pair(sz, p2))
        Tgamma(i, p2i(sz, p2)) = p2(2);
    end
    
    % exits
    Tlambda(i, i) = -p1(1);
    Tgamma(i, i)  = -p1(2);
    
end

end
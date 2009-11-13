function Pbin = bin(nn, Pb)
    Pbin = zeros(nn,1);
    for i = [1:nn]
        Pbin(i) = sum(Pb((2^(i-1)+2):(2^i+1)));
    end
end

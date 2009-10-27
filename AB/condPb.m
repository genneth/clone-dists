function Pb2 = condPb(Pb)
    Pb2 = Pb ./ (1 - Pb(0+1) - Pb(1+1));
    Pb2(0+1) = 0;
    Pb2(1+1) = 0;
end

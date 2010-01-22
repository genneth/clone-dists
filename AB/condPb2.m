function Pb2 = condPb2(Pb)
    Pb2 = Pb ./ (1 - Pb(0+1));
    Pb2(0+1) = 0;
end

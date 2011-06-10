function r = randg(A)

assert(A == int64(A));

r = sum(-log(rand([A 1])));

end
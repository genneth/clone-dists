function p = PDX(Pbs, data)

% expect Pbs to be a 2D array where Pbs(b+1,s+1) is the probability of
% finding a clone with b basal cells and s suprabasal cells.

% expect data to be a m x 3 array where each row [b s n] means clones of
% size (b,s) were seen n times

p = 1;
for d = data
    b = d(1);
    s = d(2);
    n = d(3);
    p = p * Pbs(b+1,s+1)^n;
end

end

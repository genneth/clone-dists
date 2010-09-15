function prob = binochance(s,n,p)

watermark = binopdf(s,n,p);

for low = 0:n
    if watermark < binopdf(low,n,p)
        break
    end
end
low = max(low-1,0);

for high = (low+1):n
    if watermark >= binopdf(high,n,p)
        break
    end
end

prob = 1 - (binocdf(high,n,p) - binocdf(low,n,p));

end

function [tmp] = findtrace(a, d1, d2)
%Using the topography as the trace.
n = size(a,1);
ntag = 0;
for i = 1: n 
    if a(i,3) > d1 && a(i,3) < d2
        ntag = ntag + 1;
        tmp(ntag,:) = a(i,:);
    end
    
end


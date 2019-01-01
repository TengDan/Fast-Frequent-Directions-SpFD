function [ E, V, time ] = freqDir( A, l)
% condition: l has to be at most half of d, otherwise the econ svd below
% will have problem
tic;
[n, d] = size(A);
B = A(1:l,:);
B = [B; zeros(l,d)];
cat = ceil(n/l-1);

for i = 1: cat
    if i == cat
        rp = n - (i*l) + l;
        B((l+1): rp, :) = A(i*l+1: n, :);
    else
        B((l+1): (2*l), :) = A(i*l+1: (i+1)*l, :);
    end
    [~, E, V] = svd(B, 'econ');
    E = sqrt(max(E.^2-eye(2*l)*(E(l+1, l+1)^2), 0));
    B = E * V';
end
time = toc;
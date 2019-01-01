function [ E, V, time ] = randFreqDirP( A, l, r )
tic;
[n, d] = size(A);
A = A(randperm(n), :);
B = Sparse(A(1:r,:), l);
B = [B; zeros(l, d)];
cat = ceil(n/r-1);
if cat == 0
    [~, E, V] = svd(B, 'econ');
    E = sqrt(max(E.^2-eye(2*l)*(E(l+1, l+1)^2), 0));
    B = E * V';
end
for i = 1: cat
    if i == cat
        B((l+1): (2*l), :) = Sparse(A(i*r+1: n, :), l);
    else
        B((l+1): (2*l), :) = Sparse(A(i*r+1: (i+1)*r, :), l);
    end
    [~, E, V] = svd(B, 'econ');
    E = sqrt(max(E.^2-eye(2*l)*(E(l+1, l+1)^2), 0));
    B = E * V';
end
time = toc;

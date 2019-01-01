function [SAk, time] = lowrankV(A, V, k)
tic;
if size(V, 2)<k
    fprintf('Wrong Input, k exceeds the number of rows of V!');
else
    [u, s, v] = svd(A*V, 'econ');
    SAk = u(:,1:k) * s(1:k,1:k) * v(:,1:k)' * V';
end
time = toc;

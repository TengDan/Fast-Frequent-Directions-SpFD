function [ SAk, time1, timet ] = randlowrank( A, SA, k)
tic;
if size(SA, 1) < k 
    fprintf('Wrong Input, k exceeds the number of rows of S!');
else
    [V, ~] = qr(SA',0); 
    time1 = toc;
    [U1, S1, V1] = svd(A * V, 'econ');
    SAk = U1(:,1:k) * S1(1:k,1:k) * V1(:,1:k)' * V';
end
timet = toc;


% norm sampling based on A's row norm

function [SA, time] = normSamp(A, r)
tic;
[n, d] = size(A);
normScores = sum(A.^2, 2);
prob = normScores/sum(normScores);
s = randsample(n, r, true, prob);

SA = zeros(r, d);
for i = 1:r
    SA(i,:) = A(s(i),:)/sqrt(r*prob(s(i)));
end
time = toc;
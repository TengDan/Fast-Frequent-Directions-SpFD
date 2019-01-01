% Subsampled Randomized Hadamard Transform
% From C: n x (d+m), obtain the new smaller sized D: s x (d+m)

function [D, time] = DisCT(C,s)        % Here C is a nx(d+m) matrix, [A B].
tic;
n = size(C,1);
sgn = randi(2,[n,1])*2-3;     % generate +-1 column vector
C = bsxfun(@times, C, sgn);   % elementwise mult, each row of A mult +-1

%n = 2^(ceil(log2(n)));        % select n to the smallest power of 2.
D = dct(C, n);
idx = sort(randsample(n, s));   % randomly select the rows of A
D = D(idx,:);
D = D*(n/sqrt(s));            % sqrt(n/s),fwht: /n, not/sqrt(n),compensate
                              % the n here should be a power of 2, not the
                              % original size.
time = toc;
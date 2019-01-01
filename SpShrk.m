function BB = SpShrk(rows, cols, values, n, d, l)
G = randn(d, l);
nnz = size(rows,1);
C = zeros(n,l);
for j = 1:1:nnz
    C(rows(j), :) = values(j) * G(cols(j), :) + C(rows(j), :);            
end

D = zeros(d,l);
for j = 1:1:nnz
    D(cols(j), :) = values(j) * C(rows(j), :) + D(cols(j), :); 
end
C = zeros(n,l);
for j = 1:1:nnz
    C(rows(j), :) = values(j) * D(cols(j), :) + C(rows(j), :); 
end
[Z, ~] = qr(C,0); % Z is nxl 
D = zeros(d, l);
for j = 1:1:nnz
    D(cols(j),:) = values(j)*Z(rows(j),:) + D(cols(j),:);
end

[~, E, V] = svd(D', 'econ');
BB = E*V';

% This algorithm is not suitable for dense matrix as n will be smaller than
% l. This happens when the number of nonzeros is too large. Then for each
% block, n is much smaller than d, then much smaller than l. 

% [rows, cols, values] = find(AA');

%nnz = size(rows,1);
%P = zeros(d,l);
%for j = 1:1:nnz
%    P(rows(j), :) = values(j) * Z(cols(j), :) + P(rows(j), :);            
%end
%P = P';
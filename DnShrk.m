function [B, V] = DnShrk(A, l)
[~, E, V] = svd(A, 'econ');
if size(E, 1) > l
    E = sqrt(max(E.^2-eye(size(E))*E(l+1, l+1)^2, 0));
    B = E(1:l,1:l)*V(:,1:l)';
else
    B = A;
end
